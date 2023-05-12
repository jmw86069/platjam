
#' Import Nanostring csv file
#'
#' Import Nanostring csv file
#'
#' @family jam import functions
#'
#' @return `SummarizedExperiment`
#'
#' @param csv `character` path to CSV file, suitable for import
#'    by `data.table::fread()`, so most delimiters will be recognized.
#' @param probe_colname `character` string used to define the probe
#'    identifier in `rowData()`, derived from the `colnames()` of the
#'    import data. Sometimes the values represent `"Probe_ID"`, and
#'    sometimes they may represent `"Gene"`, `"miRNA"` or another
#'    measurement.
#' @param probe_anno_file `character` path to one or more files with
#'    additional annotation for each measurement, with identifiers that
#'    match `probe_colname` above.
#' @param assay_name `character` string of optional assay name used
#'    in the `SummarizedExperiment` object created. When `NULL` the
#'    default is `"raw"` unless the file contains substring `"norm"`
#'    in which case it uses `"norm"`.
#' @param hk_count `integer` number of housekeeper genes to use when
#'    `"control_type"` is not defined for the imported data.
#'    In this case, the last `hk_count` genes in the data are
#'    assumed to be housekeeper genes, by typical convention of
#'    NanoString codeset design.
#' @param curation_txt either `data.frame` or `character` file path
#'    to tab- or comma-delimited file. The first column should match the
#'    column headers after importing data, `colData(se)`.
#'    Subsequent columns contain associated sample annotations.
#'    For Nanostring data, the Nanostring sample annotations will already
#'    be associated with the `colData(se)`, and `colnames(curation_df)`
#'    will overwrite any that already exist.
#'    Pro tip: The first column in `curation_txt` should contain `'.'`
#'    instead of punctuation/whitespace, to improve pattern matching
#'    filenames where the punctuation characters may have been modified
#'    during processing.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
import_nanostring_csv <- function
(csv,
 probe_colname="Probe_ID",
 probe_anno_file=NULL,
 assay_name=NULL,
 hk_count=10,
 curation_txt=NULL,
 verbose=FALSE,
 ...)
{
   # read csv into data.frame
   df <- data.table::fread(file=csv,
      data.table=FALSE,
      ...)
   rownames(df) <- df[,1];

   # detect sample colnames
   sample_col_idx <- grep("^(Sample|Lane)_",
      ignore.case=TRUE,
      colnames(df))
   sample_colnames <- colnames(df)[seq_len(max(sample_col_idx))];
   probe_names <- setdiff(colnames(df), sample_colnames);
   probe_df <- data.frame(check.names=FALSE,
      Probe_ID=probe_names)
   colnames(probe_df) <- probe_colname;

   # optional probe annotation file
   if (length(probe_anno_file) > 0 && all(file.exists(probe_anno_file))) {
      for (ifile in probe_anno_file) {
         idf <- data.table::fread(ifile, data.table=FALSE)
         idf_colnames <- setdiff(colnames(idf), colnames(probe_df));
         probe_df[,idf_colnames] <- idf[match(probe_names, idf[,1]), idf_colnames,drop=FALSE];
      }
   }

   # NanoString expression
   # default is to transform with log(1 + x)
   m <- log2(1 + t(as.matrix(df[, probe_names, drop=FALSE])))

   # SummarizedExperiment
   nano_se <- SummarizedExperiment::SummarizedExperiment(
      assays=list(raw=m),
      rowData=probe_df,
      colData=df[,sample_colnames, drop=FALSE])

   # optional curation_txt
   if (length(curation_txt) > 0) {
      # if curation_txt is supplied, use it to annotate the samples
      #
      if (verbose) {
         jamba::printDebug("import_nanostring_csv(): ",
            "processing curation_txt.");
      }
      if (is.atomic(curation_txt)) {
         if (length(curation_txt) == 1) {
            curation_txt <- data.table::fread(file=curation_txt,
               data.table=FALSE);
         } else {
            stop("curation_txt must be a single file path, or a data.frame")
         }
      }
      if ("data.frame" %in% class(curation_txt)) {
         # consider order_priority="x" as an option
         new_sample_df <- curate_to_df_by_pattern(
            x=colnames(nano_se),
            input_colname=head(colnames(curation_txt), 1),
            df=curation_txt,
            verbose=verbose);
         if (ncol(new_sample_df) > 1) {
            sample_match <- match(colnames(nano_se), new_sample_df$File)
            if (all(is.na(sample_match))) {
               jamba::printDebug("import_nanostring_csv(): ",
                  c("No values matched curation_txt, ",
                  "try replace non-alphanumeric characters with "), "'.'")
            }
            add_colnames <- colnames(new_sample_df);
            if (verbose) {
               jamba::printDebug("import_nanostring_csv(): ",
                  "added columns: ", add_colnames);
               # jamba::printDebug("new_sample_df:");
               # print(new_sample_df);
            }
            colData(nano_se)[,add_colnames] <- (
               new_sample_df[sample_match, add_colnames, drop=FALSE])
         }
      } else {
         stop("curation_txt must be a single file path, or a data.frame")
      }
   }

   if (length(assay_name) > 0) {
      names(SummarizedExperiment::assays(nano_se)) <- head(assay_name, 1);
   } else if (any(grepl("norm", ignore.case=TRUE, csv))) {
      names(SummarizedExperiment::assays(nano_se)) <- "norm";
   }

   # add probe control_type
   if (!"control_type" %in% colnames(SummarizedExperiment::rowData(nano_se))) {
      SummarizedExperiment::rowData(nano_se)$control_type <- ifelse(
         grepl("^POS", rownames(nano_se)),
         "POS",
         ifelse(
            grepl("^NEG", rownames(nano_se)),
            "NEG",
            NA))
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
         levels=intersect(c("NEG", "POS", "HK"),
            SummarizedExperiment::rowData(nano_se)$control_type))
      control_probes <- split(rownames(nano_se),
         SummarizedExperiment::rowData(nano_se)$control_type)
   }

   # update colnames(nano_se) to remove whitespace and file extension
   colnames(nano_se) <- gsub("[ _]+", "_",
      gsub("[.]RCC$", ignore.case=TRUE, "", colnames(nano_se)))

   return(nano_se)
}
