
#' Import Nanostring RCC files
#'
#' Import Nanostring RCC files
#'
#' This function reads one or more Nanostring RCC files,
#' and produces a `SummarizedExperiment` object.
#'
#' This function was based upon the function
#' `NanoStringNorm::read.markup.RCC()`, which in our experience
#' was problematic because it assumed all files contained identical
#' sample annotation content.
#'
#' For examples how to curate sample names into proper `data.frame`
#' of experimental factors, see `splicejam::curateVtoDF()` or
#' `splicejam::curateDFtoDF()`.
#'
#' @family jam import functions
#'
#' @return `SummarizedExperiment` or `NanoString` object, based
#' upon argument `return_type`.
#'
#' @param rcc_files character vector of file paths to RCC files,
#'    or if `NULL` then the `rcc_path` and `rcc_pattern` arguments
#'    are used to find RCC files.
#' @param rcc_path character vector of file directories to look
#'    for RCC files, the files are matched to pattern `rcc_pattern`.
#'    Note the `rcc_path` is only used when `rcc_files` is not supplied.
#' @param rcc_pattern character string containing regular expression
#'    pattern, used to match filenames which should be considered RCC
#'    files.
#' @param exclude optional character vector of filenames to exclude from
#'    import. The filenames may be either a full file path, or
#'    the basename without the file path.
#' @param include optional character vector of filenames to include
#'    during data import. Note that when `include` is defined, *only*
#'    these files are included in the file import, all other files are
#'    ignored.
#' @param return_type `character` string indicating the format to return.
#'    SummarizedExperiment is recommended, as the NanoString format
#'    was useful for R analysis in a package that is no longer
#'    actively maintained.
#' @param hk_count `integer` number of housekeeper genes to use when
#'    `"control_type"` is not defined for the imported data.
#'    In this case, the last `hk_count` genes in the data are
#'    assumed to be housekeeper genes, by typical convention of
#'    NanoString codeset design.
#' @param assay_name `character` string of optional assay name used
#'    in the `SummarizedExperiment` object created. When `NULL` the
#'    default is `"raw"`.
#' @param curation_txt `data.frame` whose first column should match the
#'    sample column headers found in the PD abundance columns, and
#'    subsequent columns contain associated sample annotations.
#'    If `curation_txt` is not supplied, then values will be split into
#'    columns by `_` underscore or `" "` whitespace characters.
#' @param debug logical indicating whether to send intermediate data
#'    before full processing, useful for debugging file format errors.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
import_nanostring_rcc <- function
(rcc_files=NULL,
 rcc_path=".",
 rcc_pattern=NULL,
 exclude=NULL,
 include=NULL,
 nprobes=-1,
 control_greps=c(POS="^POS", NEG="^NEG"),
 hk_count=10,
 assay_name=NULL,
 curation_txt=NULL,
 return_type=c("SummarizedExperiment", "NanoString"),
 debug=FALSE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to correct a bug in the importer which assumed
   ## all RCC files use the same order of sample annotations,
   ## otherwise it appends data without regard to ordering the
   ## columns consistently.
   ##
   ## If rcc_files supplies a vector of file names including proper file path,
   ## they will be used as-is. If this vector is empty, then rcc_path will
   ## be used to supply one or more directories inside which RCC files
   ## matching rcc_pattern will be used.
   ##
   ## curateFrom and curateTo are used by curateText() which runs
   ## gsub(from, to) on each line of XML read from the RCC file, prior
   ## to parsing the XML. This step was implemented to allow correcting
   ## invalid XML characters.
   ##
   return_type <- match.arg(return_type);
   if ("SummarizedExperiment" %in% return_type &&
         !suppressPackageStartupMessages(require(SummarizedExperiment))) {
      stop("The SummarizedExperiment package is required for SummarizedExperiment output.");
   }

   # define suitable default rcc_pattern if needed
   if (is.null(rcc_pattern)) {
      rcc_pattern <- "[.](RCC|rcc)$";
   }
   if (!suppressPackageStartupMessages(require(XML))) {
      stop("The XML package is required for import_nanostring_rcc() to parse the XML result files.");
   }


   if (length(rcc_files) > 0 && any(file.exists(rcc_files))) {
      if (any(!file.exists(rcc_files))) {
         jamba::printDebug("import_nanostring_rcc(): ",
            "Error: Some (",
            sum(!file.exists(rcc_files)),
            ") supplied rcc_files do not exist. Please remedy.");
         stop("rcc_files not found.");
      }
   } else {
      if (!file.exists(rcc_path)) {
         stop("rcc_path must be a valid file path.") ;
      }
      if (verbose) {
         jamba::printDebug("import_nanostring_rcc(): ",
            "rcc_pattern:",
            rcc_pattern);
      }
      rcc_files <- list.files(path=rcc_path,
         pattern=rcc_pattern,
         full.names=TRUE);
   }
   if (length(exclude) > 0) {
      rcc_exclude <- (rcc_files %in% exclude |
            basename(rcc_files) %in% exclude);
      if (verbose && any(rcc_exclude)) {
         jamba::printDebug("import_nanostring_rcc(): ",
            "excluding ",
            sum(rcc_exclude),
            " files matching exclude string.");
      }
      rcc_files <- rcc_files[!rcc_exclude];
   }
   if (length(include) > 0) {
      rcc_include <- (rcc_files %in% include |
            basename(rcc_files) %in% include);
      if (verbose && any(!rcc_include)) {
         jamba::printDebug("import_nanostring_rcc(): ",
            "including only ",
            sum(rcc_include),
            " files matching include string.");
      }
      rcc_files <- rcc_files[rcc_include];
   }


   # check if any rcc_files
   if (length(rcc_files) == 0) {
      jamba::printDebug("import_nanostring_rcc(): ",
         "No files were found in rcc_path:", head(rcc_path),
         ", rcc_pattern:", head(rcc_pattern),
         ", exclude:", exclude,
         ", include:", include);
      stop("import_nanostring_rcc(): no RCC files available.");
   }

   rcc.header.merged <- NULL;
   rcc.data.merged <- NULL;

   ## Switch to XML parsing
   #count <- 1;
   #for (rcc_file in rcc_files) {
   for (count in seq_along(rcc_files)) {
      rcc_file <- rcc_files[[count]];
      if (verbose) {
         jamba::printDebug("import_nanostring_rcc(): ",
            "reading rcc_file:",
            rcc_file);
      }
      clean_nano_filename <- function(x){
         gsub("[.]RCC$",
            "",
            ignore.case=TRUE,
            gsub(" ",
               "_",
               basename(x)));
      }
      sample_name <- clean_nano_filename(rcc_file);
      #rccXml2 <- xmlToDataFrame(c("<RCC>", readLines(rcc_file), "</RCC>"));
      #rccXml1 <- xmlTreeParse(c("<RCC>", readLines(rcc_file), "</RCC>"));
      #rccXml <- xmlRoot(rccXml1);

      ## Parse XML
      rccLines <- readLines(rcc_file);
      #rccLines <- curateText(readLines(rcc_file),
      #   curateFrom=curateFrom,
      #   curateTo=curateTo);
      rccXmlL <- xmlToList(xmlRoot(xmlTreeParse(
         c("<RCC>",
            rccLines,
            "</RCC>"))));
      if (verbose) {
         jamba::printDebug("import_nanostring_rcc(): ",
            "parsed XML");
      }

      headerNames <- jamba::provigrep(c("header","sample","lane"), names(rccXmlL));
      codeNames <- jamba::vigrep("code", names(rccXmlL));
      if (verbose) {
         jamba::printDebug("import_nanostring_rcc(): ",
            "headerNames:", headerNames);
         jamba::printDebug("import_nanostring_rcc(): ",
            "codeNames:", codeNames);
      }
      rcc.header <- jamba::rbindList(lapply(jamba::nameVector(headerNames), function(i){
         iXml <- rccXmlL[[i]];
         iDF <- data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            jamba::rbindList(
               strsplit(unlist(
                  strsplit(iXml, "[\r\n]")
               ), ",")
            )
         );
         ## For Sample_Attributes, Lane_Attributes, we create colnames
         ## Sample_*, and Lane_*
         if (jamba::igrepHas("_Attributes", i)) {
            iDF1 <- paste0(gsub("_Attributes.*$", "_", i), iDF[,1]);
            iDF[,1] <- iDF1;
         }
         colnames(iDF)[2] <- clean_nano_filename(rcc_file);
         iDF;
      }));
      if (verbose) {
         jamba::printDebug("import_nanostring_rcc(): ",
            "parsed header");
         jamba::printDebug("import_nanostring_rcc(): ",
            "codeNames:", codeNames);
      }
      if (debug) {
         return(rccXmlL);
      }
      rcc.data1 <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         jamba::rbindList(
            strsplit(unlist(
               strsplit(rccXmlL[[codeNames]], "\n")
            ), ",")
         )
      );
      rcc.data <- read.table(
         text=jamba::pasteByRow(sep="\t", rcc.data1),
         stringsAsFactors=FALSE,
         sep="\t",
         header=TRUE,
         comment.char="");
      #rcc.data <- data.frame(
      #   nameList(
      #      as.list(rcc.data1[-1,,drop=FALSE]),
      #      unlist(rcc.data1[1,])));
      rcc.data <- jamba::renameColumn(rcc.data,
         from="Count",
         to=sample_name);
      if (nprobes > 0) {
         if (verbose) {
            jamba::printDebug("import_nanostring_rcc(): ",
               "loading only the first ",
               nprobes,
               " probes.");
         }
         rcc.data <- head(rcc.data, nprobes);
      }
      if (verbose) {
         jamba::printDebug("import_nanostring_rcc(): ",
            "parsed data");
      }

      if (count == 1) {
         rcc.header.merged <- rcc.header;
         rcc.data.merged <- rcc.data;
      } else {
         rcc.header.merged <- data.frame(rcc.header.merged,
            subset(rcc.header[match(rcc.header.merged[,1], rcc.header[,1]),,drop=FALSE], select=2),
            check.names=FALSE);
         rcc.data.merged <- data.frame(rcc.data.merged, subset(rcc.data, select=4),
            check.names=FALSE);
      }

      #count <- count + 1;
   }

   #rcc.header.merged[3,1] <- "sample.id";
   #rcc.header.merged[9,1] <- "lane.id";
   rownames(rcc.header.merged) <- jamba::makeNames(rcc.header.merged[,1]);
   rcc.header.merged <- rcc.header.merged[,-1,drop=FALSE];

   #colnames(rcc.data.merged[1]) <- "Code.Count";

   x <- list(x=rcc.data.merged,
      header=rcc.header.merged);

   if ("NanoString" %in% return_type) {
      class(x) <- 'NanoString';
      return(x);
   }

   ## Convert to SummarizedExperiment
   rownames(x$x) <- jamba::makeNames(x$x$Name);
   gene_colnames <- setdiff(colnames(x$x), colnames(x$header));
   nano_assays <- as.matrix(x$x[,colnames(x$header),drop=FALSE]);

   # by default transform with log2(1 + x)
   nano_assays <- log2(1 + nano_assays);

   nano_genes <- x$x[,gene_colnames,drop=FALSE];
   if (length(control_greps) > 0) {
      nano_controls_l <- jamba::provigrep(control_greps,
         rownames(nano_genes),
         returnType="list");
      nano_controls <- jamba::nameVector(
         rep(names(nano_controls_l), lengths(nano_controls_l)),
         unlist(nano_controls_l)
      )
      nano_genes$control_type <- NA;
      nano_genes[names(nano_controls),"control_type"] <- nano_controls;
   }

   # parse sample annotations from Nanostring
   nano_samples <- read.table(
      text=c(paste(rownames(x$header), collapse="\t"),
         jamba::pasteByRow(sep="\t",
            condenseBlanks=FALSE,
            t(x$header))),
      stringsAsFactors=FALSE,
      sep="\t",
      header=TRUE,
      comment.char="");
   rownames(nano_samples) <- colnames(x$header);

   # import sample annotations
   isamples <- rownames(nano_samples);
   if (length(curation_txt) > 0) {
      if (verbose) {
         jamba::printDebug("import_nanostring_rcc(): ",
            "applying sample annotations via curation_txt.");
      }
      # sample_df
      sample_df <- curate_to_df_by_pattern(
         isamples,
         df=curation_txt,
         ...);
      isamples1 <- rownames(sample_df);

      # match the order to nano_samples
      filename_colname <- head(jamba::vigrep("^filename$",
         colnames(sample_df)),
         1);
      nano_sample_match <- match(sample_df[[filename_colname]],
         rownames(nano_samples))
      use_nano_samples <- nano_samples[nano_sample_match, , drop=FALSE];
      colnames(use_nano_samples) <- tail(makeNames(c(
         colnames(sample_df), colnames(nano_samples))),
         ncol(nano_samples))
      sample_df[, colnames(use_nano_samples)] <- use_nano_samples[
         nano_sample_match, , drop=FALSE];

      if (!all(isamples1 == isamples)) {
         # synchronize sample order
         isamples_match <- match(sample_df[[filename_colname]],
            isamples);
         isamples_from <- isamples[isamples_match];
         isamples_to <- isamples1;
         if (verbose > 1) {
            jamba::printDebug("import_nanostring_rcc(): ",
               "renaming sample colnames to match sample_df:");
            print(data.frame(isamples_from,
               isamples_to));
         }
         # rename matrix colnames
         nano_assays <- jamba::renameColumn(nano_assays,
            from=isamples_from,
            to=isamples_to);
      }
      isamples <- isamples1;
      if (verbose) {
         if (verbose > 1) {
            jamba::printDebug("import_salmon_quant(): ",
               "sample_df:");
            print(sample_df);
         }
      }
      nano_samples <- sample_df;
   }

   nano_se <- SummarizedExperiment::SummarizedExperiment(
      assays=list(raw=nano_assays[, isamples, drop=FALSE]),
      colData=nano_samples[match(isamples, rownames(nano_samples)), , drop=FALSE],
      rowData=nano_genes);
   if (length(assay_name) > 0) {
      names(SummarizedExperiment::assays(nano_se)) <- head(assay_name, 1);
   }

   # add probe control_type if needed
   if (!"control_type" %in% colnames(rowData(nano_se))) {
      SummarizedExperiment::rowData(nano_se)$control_type <- ifelse(
         grepl("^POS_", rownames(nano_se)),
         "POS",
         ifelse(
            grepl("^NEG_", rownames(nano_se)),
            "NEG",
            NA))
   }
   # attempt to add housekeeper genes if not defined
   if (!any(grepl("^h|control", ignore.case=TRUE, SummarizedExperiment::rowData(nano_se)$control_type))) {
      control_probes <- split(rownames(nano_se),
         SummarizedExperiment::rowData(nano_se)$control_type)
      # guess HK genes
      hk_genes <- tail(
         setdiff(rownames(nano_se), unlist(control_probes)),
         hk_count)
      hkmatch <- match(hk_genes, rownames(nano_se));
      SummarizedExperiment::rowData(nano_se[hkmatch,])$control_type <- "HK";
      SummarizedExperiment::rowData(nano_se)$control_type <- factor(
         SummarizedExperiment::rowData(nano_se)$control_type,
         levels=jamba::provigrep(c("NEG", "POS", "HK", "."),
            unique(SummarizedExperiment::rowData(nano_se)$control_type)))
      control_probes <- split(rownames(nano_se),
         SummarizedExperiment::rowData(nano_se)$control_type)
   }

   return(nano_se);
}

