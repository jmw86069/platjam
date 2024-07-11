
#' Import generic omics tsv or csv file
#'
#' Import generic omics tab- or comma-delimited data file
#'
#' @family jam import functions
#'
#' @return `SummarizedExperiment`
#'
#' @param x one of the following input types:
#'    * `character` path to a data file suitable for import
#'    by `data.table::fread()`, where most delimiters will be recognized
#'    automatically.
#'    * `data.frame` representing the equivalent data as from a file.
#'    * `matrix` representing only the measured values, with rownames
#'    and colnames that represent respective identifiers.
#' @param assay_name `character` string which will define the assay name
#'    where data is stored in the resulting `SummarizedExperiment` object.
#' @param row_identifier defines the column, columns, or rownames, to be
#'    used as row identifiers. These values become `rownames(se)` for
#'    the output object.
#'    One of the following:
#'    * `NULL` (default) which auto-detects an appropriate row identifier:
#'
#'       * If the first column is numeric, it uses rownames.
#'       * Otherwise the first column is used. If there are duplicated values,
#'       they are made unique with `jamba::makeNames()`.
#'
#'    * `character` string of colname in the data `x` imported.
#'    * `integer` column number of data `x` imported.
#'    * The integer value `0` or `-1` to indicate `rownames(x)`
#'    * One of these strings, to indicate `rownames(x)`:
#'    `"rownames"`, `"rowname"`, or `"row.names"`.
#'    * The default is to use the first column in the imported data `x`.
#' @param row_annotation_columns defines columns to retain in `rowData()`
#'    which have non-measurement data associated with each row.
#'    One of the following:
#'    * `character` vector with one or more `colnames()`
#'    * `integer` vector with one or more column numbers
#' @param curation_txt either `data.frame` or `character` file path
#'    to tab- or comma-delimited file.
#'    The first column should match the
#'    column headers after importing data, `colData(se)`.
#'    Subsequent columns contain associated sample annotations.
#'    For Nanostring data, the Nanostring sample annotations will already
#'    be associated with the `colData(se)`, and `colnames(curation_df)`
#'    will overwrite any that already exist.
#'    * Pro tip: The first column in `curation_txt` should contain `'.'`
#'    instead of punctuation/whitespace, to improve pattern matching
#'    filenames where the punctuation characters may have been modified
#'    during processing.
#'    * Note also that when `curation_txt` is supplied, samples in `se`
#'    will be subset to include only those samples that match `curation_txt`,
#'    and in the order they appear in the `curation_txt` file.
#'    This behavior allows the `curation_txt` to be used to define
#'    the appropriate experimental ordering, which by default also
#'    defines downstream control factor levels for statistical contrasts.
#'    The first factor level is used as the control value in those
#'    contrasts.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to internal functions
#'    for example `data.table::fread()`.
#'
#' @examples
#' x <- matrix(letters[1:9], ncol=3);
#' rownames(x) <- LETTERS[1:3];
#' colnames(x) <- letters[1:3];
#' x
#' import_omics_data(x)
#'
#' x <- matrix(1:9, ncol=3);
#' rownames(x) <- LETTERS[1:3];
#' colnames(x) <- paste0(letters[1:3], "_rep1");
#' x
#'
#' curation_txt <- data.frame(Pattern=LETTERS[1:3], Group="A", Batch="B")
#' se <- import_omics_data(x, curation_txt=curation_txt)
#' SummarizedExperiment::assays(se)[[1]]
#' SummarizedExperiment::colData(se)
#' SummarizedExperiment::rowData(se)
#'
#' se2 <- import_omics_data(x, curation_txt=curation_txt[c(2, 1, 3), ])
#' SummarizedExperiment::assays(se2)[[1]]
#' SummarizedExperiment::colData(se2)
#'
#' @export
import_omics_data <- function
(x,
 assay_name="data",
 row_identifier=NULL,
 row_annotation_columns=NULL,
 curation_txt=NULL,
 verbose=FALSE,
 ...)
{
   #
   if ("character" %in% class(x)) {
      # file import
      if (length(x) != 1) {
         cli::cli_abort(message=paste0(
            "{.field x} when supplied as {.cls character} must have ",
            "only one value."))
      }
      # read csv into data.frame
      df <- data.table::fread(file=x,
         data.table=FALSE,
         ...)
   } else if ("data.frame" %in% class(x)) {
      df <- x;
      # use as-is
   } else if ("matrix" %in% class(x)) {
      if (!is.numeric(x)) {
         cli::cli_abort(message=paste0(
            "When matrix data is supplied, it must be {.cls numeric} ",
            "as determined by {.code is.numeric(x)}"));
      }
      row_identifier <- -1;
      df <- x;
   } else if ("tibble" %in% class(x)) {
      # coerce to data.frame
      df <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         x);
   }

   ## check dimensions
   if (nrow(df) == 0 || ncol(df) == 0) {
      cli::cli_abort(message=paste0(
         "Imported data must have at least one column and one row, ",
         "nrow={nrow(df)}, ncol={ncol(df)}"))
   }

   ## handle empty row_identifier
   if (length(row_identifier) == 0) {
      # auto-detect
      if (is.numeric(df[[1]])) {
         row_identifier <- -1;
      } else {
         row_identifier <- 1;
      }
   }

   ## handle character row_identifier, convert to numeric
   if (is.factor(row_identifier)) {
      row_identifier <- as.character(row_identifier);
   }
   if (is.character(row_identifier)) {
      matchcol <- match(row_identifier,
         c(colnames(df), "rownames", "rowname", "row.names"));
      if (any(is.na(matchcol))) {
         cli::cli_abort(message=paste0(
            "{.var row_identifier} as character ({.val {row_identifier}}",
            ") value did not ",
            "match imported data colnames: ",
            "{.val {colnames(df)}}"))
         stop(paste0("row_identifier did not match colnames of imported data."))
      }
      if (any(matchcol > ncol(df))) {
         matchcol[matchcol > ncol(df)] <- -1;
      }
      row_identifier <- unique(matchcol);
   }
   if (!is.numeric(row_identifier)) {
      accepted_classes <- c("integer",
         "numeric",
         "character",
         "factor",
         "NULL");
      cli::cli_abort(message=paste0(
         "{.var row_identifier} accepts these classes: ",
         "{.val {accepted_classes}}"))
   }
   row_identifier[row_identifier <= 0] <- -1;
   row_identifier <- unique(row_identifier);
   row_identifiers <- jamba::makeNames(
      jamba::pasteByRow(do.call(cbind,
         lapply(row_identifier, function(ir) {
            if (ir < 0) {
               rownames(df)
            } else {
               df[[ir]]
            }
         })),
         ...),
      ...)
   rownames(df) <- row_identifiers;

   ## handle row_annotation_columns
   data_colnums <- seq_len(ncol(df));
   if (any(row_identifier > 0)) {
      data_colnums <- data_colnums[-row_identifier[row_identifier > 0]];
   }
   row_annotation_columns <- jamba::rmNA(row_annotation_columns);
   if (length(row_annotation_columns) > 0) {
      if (is.numeric(row_annotation_columns)) {
         row_annotation_columns <- unique(row_annotation_columns[
            row_annotation_columns > 0 &
            row_annotation_columns <= ncol(df)]);
         row_annotation_columns <- colnames(df)[row_annotation_columns];
      }
   }
   row_annotation_colnums <- NULL;
   if (length(row_annotation_columns) > 0) {
      row_annotation_columns <- as.character(row_annotation_columns);
      use_row_annotation_columns <- intersect(row_annotation_columns,
         colnames(df));
      if (length(use_row_annotation_columns) == 0) {
         cli::cli_abort(message=paste0(
            "{.var row_annotation_columns} ({.val {row_annotation_columns}}",
            ") did not match any data colnames: ",
            "{.val {colnames(df)}}"))
      }
      row_annotation_colnums <- match(use_row_annotation_columns,
         colnames(df));
      if (any(row_identifiers > 0)) {
         row_annotation_colnums <- unique(c(row_annotation_colnums,
            row_identifiers[row_identifiers > 0]));
      }
      data_colnums <- setdiff(data_colnums, row_annotation_colnums);
   }

   ## Verify remaining data_colnums are numeric
   data_colnums_isnum <- sapply(data_colnums, function(ic){
      is.numeric(df[[ic]])
   })
   if (any(data_colnums_isnum %in% FALSE)) {
      # add to row_annotation_colnums
      if (verbose) {
         jamba::printDebug("import_omics_data(): ",
            "data_colnums moved to row_annotation: ",
            data_colnums[!data_colnums_isnum])
      }
      row_annotation_colnums <- unique(c(row_annotation_colnums,
         data_colnums[!data_colnums_isnum]));
      data_colnums <- data_colnums[data_colnums_isnum];
   }

   ## Assemble assay matrix
   assay_matrix <- as.matrix(df[, data_colnums, drop=FALSE])

   ## Assemble rowData
   if (length(row_annotation_colnums) >= 1) {
      rowData_df <- df[, row_annotation_colnums, drop=FALSE];
   } else {
      rowData_df <- data.frame(check.names=FALSE,
         rows=row_identifiers);
      rownames(rowData_df) <- row_identifiers;
   }

   ## Prepare colData
   if (length(curation_txt) == 0) {
      colData_df <- data.frame(check.names=FALSE,
         columns=colnames(assay_matrix));
      rownames(colData_df) <- colnames(assay_matrix);
   } else {
      if (verbose) {
         jamba::printDebug("import_omics_data(): ",
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
      if (any(c("tibble", "matrix") %in% class(curation_txt))) {
         curation_txt <- data.frame(check.names=FALSE, curation_txt);
      }
      if ("data.frame" %in% class(curation_txt)) {
         # consider order_priority="x" as an option
         input_colname <- head(colnames(curation_txt), 1);
         new_sample_df <- curate_to_df_by_pattern(
            x=colnames(assay_matrix),
            input_colname=input_colname,
            df=curation_txt,
            verbose=verbose);
         if (ncol(new_sample_df) > 1) {
            # jamba::printDebug("new_sample_df:");print(new_sample_df);# debug
            # re-order columns to match curation_txt
            sample_match <- match(new_sample_df[[input_colname]],
               colnames(assay_matrix))
            # jamba::printDebug("sample_match: ", sample_match);# debug
            # Re-order, which omits samples not defined in curation_txt
            assay_matrix <- assay_matrix[, sample_match, drop=FALSE];
            colData_df <- new_sample_df;
            colnames(assay_matrix) <- rownames(colData_df);
         }
      } else {
         cli::cli_abort(message=paste0(
            "{.var curation_txt} was not usable as a ",
            "{.cls data.frame}."));
         stop("curation_txt was not usable as a data.frame")
      }
   }

   ## Create SummarizedExperiment data
   assays_list <- list(assay_matrix);
   names(assays_list) <- head(assay_name, 1);
   se <- SummarizedExperiment::SummarizedExperiment(
      assays=assays_list,
      rowData=rowData_df,
      colData=colData_df)
   return(se);
}
