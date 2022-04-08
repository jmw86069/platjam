
#' Import LipoType lipidomics data from CSV file
#'
#' Import LipoType lipidomics data from CSV file
#'
#' @family jam import functions
#'
#' @return `SummarizedExperiment` where:
#'   * `rowData` represents lipid measurement annotations
#'   * `colData` represents sample annotations, optionally including
#'   annotations via a `data.frame` supplied as `curation_txt`.
#'
#' @param csv `character` path to a comma-delimited file (csv) as
#'    provided in a LipoType summary report.
#' @param curation_txt `data.frame` whose first column should match the
#'    column headers found in the LipoType measurement data.
#'    If `curation_txt` is not supplied, then values will be split into
#'    columns by `_` underscore or `" "` whitespace characters.
#' @param ... additional arguments are ignored.
#'
#' @export
import_lipotype_csv <- function
(csv,
 curation_txt=NULL,
 num_header_rows=3,
 num_anno_columns=7,
 ...)
{
   # check csv file exists
   if (!file.exists(csv)) {
      stop(paste("csv file is not accessible:", csv));
   }

   # convenient adjustments
   num_header_rows1 <- num_header_rows - 1;
   num_anno_columns1 <- num_anno_columns + 1;

   # load sample annotations
   df1h <- data.table::fread(csv,
      skip=0,
      nrows=num_header_rows1,
      data.table=FALSE,
      header=TRUE);
   isamples <- tail(colnames(df1h), -num_anno_columns1);
   colData_csv <- data.frame(check.names=FALSE,
      t(df1h[,isamples,drop=FALSE]));
   colnames(colData_csv) <- t(df1h)[num_anno_columns1,,drop=TRUE];
   id_rowname <- colnames(df1h)[num_anno_columns1];
   id_colnames <- c(id_rowname,
      colnames(colData_csv));
   colData_csv[[id_rowname]] <- rownames(colData_csv);
   colData_csv <- colData_csv[, id_colnames, drop=FALSE];

   # load row annotations and measurements
   df1 <- data.table::fread(csv,
      skip=num_header_rows,
      data.table=FALSE,
      header=FALSE);
   rownames(df1) <- df1[[1]];
   colnames(df1) <- colnames(df1h);

   rowData_csv <- df1[,seq_len(num_anno_columns), drop=FALSE];

   # generate numeric matrix for measurement data
   imatrix <- as.matrix(df1[,isamples, drop=FALSE]);

   # create SummarizedExperiment object
   LipidSE <- SummarizedExperiment::SummarizedExperiment(
      assays=list(raw=imatrix),
      rowData=rowData_csv,
      colData=colData_csv)
   return(LipidSE);
}
