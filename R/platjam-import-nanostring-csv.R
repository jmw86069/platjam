
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
#' @param hk_count `integer` number of housekeeper genes to use when
#'    `"control_type"` is not defined for the imported data.
#'    In this case, the last `hk_count` genes in the data are
#'    assumed to be housekeeper genes, by typical convention of
#'    NanoString codeset design.
#' @param ... additional arguments are ignored.
#'
#' @export
import_nanostring_csv <- function
(csv,
 probe_colname="Probe_ID",
 probe_anno_file=NULL,
 assay_name=NULL,
 hk_count=10,
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

   # add probe control_type
   if (!"control_type" %in% colnames(rowData(nano_se))) {
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
