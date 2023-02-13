
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
 control_greps=c(POS="^POS_", NEG="^NEG_"),
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
         printDebug("import_nanostring_rcc(): ",
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
         printDebug("import_nanostring_rcc(): ",
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
         printDebug("import_nanostring_rcc(): ",
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
         printDebug("import_nanostring_rcc(): ",
            "including only ",
            sum(rcc_include),
            " files matching include string.");
      }
      rcc_files <- rcc_files[rcc_include];
   }


   # check if any rcc_files
   if (length(rcc_files) == 0) {
      printDebug("import_nanostring_rcc(): ",
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
         printDebug("import_nanostring_rcc(): ",
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
         printDebug("import_nanostring_rcc(): ",
            "parsed XML");
      }

      headerNames <- jamba::provigrep(c("header","sample","lane"), names(rccXmlL));
      codeNames <- jamba::vigrep("code", names(rccXmlL));
      if (verbose) {
         printDebug("import_nanostring_rcc(): ",
            "headerNames:", headerNames);
         printDebug("import_nanostring_rcc(): ",
            "codeNames:", codeNames);
      }
      rcc.header <- jamba::rbindList(lapply(nameVector(headerNames), function(i){
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
         printDebug("import_nanostring_rcc(): ",
            "parsed header");
         printDebug("import_nanostring_rcc(): ",
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
         text=pasteByRow(sep="\t", rcc.data1),
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
            printDebug("import_nanostring_rcc(): ",
               "loading only the first ",
               nprobes,
               " probes.");
         }
         rcc.data <- head(rcc.data, nprobes);
      }
      if (verbose) {
         printDebug("import_nanostring_rcc(): ",
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

   nano_samples <- read.table(
      text=c(paste(rownames(x$header), collapse="\t"),
         pasteByRow(sep="\t",
            condenseBlanks=FALSE,
            t(x$header))),
      stringsAsFactors=FALSE,
      sep="\t",
      header=TRUE,
      comment.char="");
   rownames(nano_samples) <- colnames(x$header);
   nano_se <- SummarizedExperiment::SummarizedExperiment(
      assays=list(raw=nano_assays),
      colData=nano_samples,
      rowData=nano_genes);
   if (length(assay_name) > 0) {
      names(SummarizedExperiment::assays(nano_se)) <- head(assay_name, 1);
   }
   return(nano_se);
}

