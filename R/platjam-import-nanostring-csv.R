
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
#' @param ... additional arguments are ignored.
#'
#' @export
import_nanostring_csv <- function
(csv,
 probe_colname="Probe_ID",
 probe_anno_file=NULL,
 assay_name=NULL,
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

   if (length(assay_name) > 0) {
      names(SummarizedExperiment::assays(nano_se)) <- head(assay_name, 1);
   } else if (any(grepl("norm", ignore.case=TRUE, csv))) {
      names(SummarizedExperiment::assays(nano_se)) <- "norm";
   }

   return(nano_se)
}
