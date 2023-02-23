
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
 type="raw",
 curation_txt=NULL,
 num_header_rows=3,
 num_anno_columns=7,
 verbose=FALSE,
 ...)
{
   if (length(csv) > 1) {
      SE_list <- lapply(seq_along(csv), function(i) {
         icsv <- csv[i];
         if (length(names(icsv)) > 0) {
            itype <- names(icsv);
         } else if (grepl("molp", icsv, ignore.case=TRUE)) {
            itype <- "molp";
         } else if (grepl("pmol", icsv, ignore.case=TRUE)) {
            itype <- "pmol";
         } else {
            itype <- icsv
         }
         if (verbose) {
            jamba::printDebug("import_lipotype_csv(): ",
               "icsv: ", icsv);
         }
         import_lipotype_csv(csv=icsv,
            type=itype,
            curation_txt=curation_txt,
            num_header_rows=num_header_rows,
            num_anno_columns=num_anno_columns,
            verbose=verbose,
            ...);
      });
      SE_names <- sapply(SE_list, function(iSE){
         names(assays(iSE))
      });
      names(SE_list) <- SE_names;
      identical_rowData <- sapply(tail(SE_list, -1), function(iSE){
         identical(
            SummarizedExperiment::rowData(SE_list[[1]]),
            SummarizedExperiment::rowData(iSE))
      })
      identical_colData <- sapply(tail(SE_list, -1), function(iSE){
         identical(
            SummarizedExperiment::colData(SE_list[[1]]),
            SummarizedExperiment::colData(iSE))
      });
      if (all(identical_rowData) && all(identical_colData)) {
         LipoTypeSE <- SE_list[[1]];
         for (iname in tail(SE_names, -1)) {
            SummarizedExperiment::assays(LipoTypeSE)[[iname]] <- SummarizedExperiment::assays(SE_list[[iname]])[[iname]];
         }
         return(LipoTypeSE)
      }
      jamba::printDebug("import_lipotype_csv(): ",
         "Something went wrong during import of multiple csv files. rowData or colData did not match.");
      return(SE_list);
   }
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

   # optionally load curation_txt
   if (length(curation_txt) > 0) {
      if (verbose) {
         jamba::printDebug("import_lipotype_csv(): ",
            "applying sample annotations via curation_txt.");
      }
      sample_df <- slicejam::curate_to_df_by_pattern(
         isamples,
         df=curation_txt,
         ...);
      filename_column <- vigrep("filename", colnames(sample_df));
      if (length(filename_column) == 0) {
         filename_column <- tail(colnames(sample_df), 1)
      }
      curation_match <- match(sample_df[[filename_column]],
         rownames(colData_csv));
      if (verbose > 1) {
         jamba::printDebug("import_lipotype_csv(): ",
            "sample_df:");
         print(sample_df);
         jamba::printDebug("import_lipotype_csv(): ",
            "curation_match: ", curation_match);
      }
      colData_csv <- colData_csv[curation_match,,drop=FALSE];
      isamples <- jamba::nameVector(rownames(colData_csv),
         rownames(sample_df));
      if (verbose > 1) {
         jamba::printDebug("import_lipotype_csv(): ",
            "isamples:");
         print(data.frame(isamples));
      }
      colData_csv[,colnames(sample_df)] <- sample_df;
   }


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

   # check for renaming isamples
   if (length(names(isamples)) > 0 &&
         !all(isamples == names(isamples))) {
      new_rownames <- names(isamples)[match(rownames(colData_csv), isamples)];
      if (verbose > 1) {
         jamba::printDebug("import_lipotype_csv(): ",
            "Renaming colnames(imatrix) and rownames(colData)");
         print(data.frame(isamples, new_rownames));
      }
      imatrix <- jamba::renameColumn(imatrix,
         from=isamples,
         to=names(isamples));
      if (any(is.na(new_rownames))) {
         jamba::printDebug("import_lipotype_csv(): ",
            fgText=c("darkorange2", "firebrick3"),
            "Error during isamples match with rownames(colData), usually caused by curation_txt:");
         stop("rownames(colData) mismatch with isamples.")
      }
      rownames(colData_csv) <- new_rownames;
      isamples <- rownames(colData_csv);
   }

   # create SummarizedExperiment object
   assay_list <- list(imatrix[,isamples, drop=FALSE]);
   names(assay_list) <- type;
   LipidSE <- SummarizedExperiment::SummarizedExperiment(
      assays=assay_list,
      rowData=rowData_csv,
      colData=colData_csv[isamples, , drop=FALSE]);
   return(LipidSE);
}
