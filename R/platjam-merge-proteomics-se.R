
#' merge proteomics SE objects
#'
#' merge proteomics SE objects
#'
#' See notes for specific arguments for a description of how
#' data is merged relative to
#' rows and `rowData()`, columns and `colData()`.
#'
#' The general strategy is to merge equivalent rows to integrate rows
#' across `SE1` and `SE2`, but to force columns (sample measurements)
#' to be unique across `SE1` and `SE2`.
#'
#' This process is somewhat similar to calling `cbind()`, in that
#' the sample columns are extended. However, the rows are merged where
#' possible.
#'
#' No assay measurement values are lost during this process.
#'
#' @family jam utility functions
#'
#' @param SE1,SE2 `SummarizedExperiment` objects to be merged into
#'    one output object.
#' @param rowname1,rowname2 `character` string that describes which
#'    `SummarizedExperiment::rowData()` annotation to use to create
#'    appropriate rownames to be merged. This approach is useful when
#'    merging data based upon gene symbol, instead of a protein accession
#'    or peptide sequence. The intent is to allow "equivalent" rows
#'    to be combined across `SE1` and `SE2`, while non-equivalent rows
#'    unique to `SE1` or `SE2` are represented on their own row.
#'
#'    The default values assume each proteomics
#'    SE object contains a rowData column `"SYMBOL"` with the official
#'    gene symbol represented on each row. This column is appropriate if
#'    proteomics data already represents abundance measurements which
#'    were already aggregated to the protein-level (i.e. gene locus level).
#'    The data will therefore be merged based upon the gene symbol.
#'    In the event that multiple rows represent the same gene symbol,
#'    they will be renamed using `jamba::makeNames(..., renameFirst=FALSE)`
#'    so that the entries will be merged in order they appear in each
#'    dataset.
#'
#'    However, if the input data contains peptide-level measurements,
#'    the appropriate column should contain the peptide sequence, so that
#'    the data is merged based upon equivalent peptide sequences.
#'
#'    If `rowname1` or `rowname2` contain multiple values, and/or are
#'    not equal to each other, a new column `"merge_key"` is created
#'    in both `SE1` and `SE2`, and populated with relevant values.
#'    When multiple columns are indicated, they are concatenated
#'    using `jamba::pasteByRow()` to fill the column `"merge_key"`.
#'    Then both `rowname1` and `rowname2` are redefined to
#'    `"merge_key"`. Note that any pre-existing `"merge_key"` column
#'    will be overwritten.
#'
#'    A combination of `"rownames"` and `colnames(rowData())` can
#'    be used.
#'
#'    The argument value should contain one value from either:
#'    1. `colnames(rowData())` for the relevant object `SE1` or `SE2`,
#'    representing a row annotation to use as the merge key. Note that
#'    any empty values (`NA` or blank string `""`) will be replaced
#'    by existing `rownames()`.
#'    2. `"rownames"` to indicate that existing `rownames()` of the
#'    relevant object `SE1` or `SE2` should be used as the merge key.
#'    Note that if a column `"rownames"` already exists in `rowData()`
#'    it will be used as-is.
#' @param rowData_colnames_intersect,colData_colnames_intersect `logical`
#'    indicating whether to retain only the intersection of
#'    `colnames(rowData())` and `colnames(colData())` in the output
#'    rowData and colData, respectively.
#'    * `TRUE`: only the intersection is retained in the output data, default.
#'    * `FALSE`: not yet implemented.
#' @param rowData_colnames_unique `character` vector with optional
#'    `colnames(rowData())` which should be retained in a uniquely-named
#'    output column, to keep its values distinct between `SE1` and `SE2`.
#'    This argument is useful for something like `"score"` where independent
#'    datasets are expected to have unique values, and which may be
#'    important to compare.
#'    Note that columns not already being retained will be ignored.
#' @param assay_names `character` vector with one or more specific assay
#'    names to retain in the output data. By default, all assay names
#'    are retained.
#' @param se_names `character` vector length=2 to define the output labels
#'    used to indicate which rows and columns were present in `SE1` and `SE2`.
#' @param startN `integer` number passed to `jamba::makeNames()` to define
#'    the suffix number for the first versioned output. Note that
#'    `renameFirst=FALSE` so the first occurrence of a `character` string
#'    will not be renamed. When `startN=2`, subsequent repeated
#'    entries will have suffix `"_v2"`, then `"_v3"` and so on.
#' @param ... additional arguments are passed to `jamba::makeNames()`.
#'
#' @export
merge_proteomics_se <- function
(SE1,
 SE2,
 rowname1="SYMBOL",
 rowname2="SYMBOL",
 rowData_colnames_intersect=TRUE,
 colData_colnames_intersect=TRUE,
 rowData_colnames_unique=c("percentCoverage",
    "numPepsUnique",
    "scoreUnique"),
 assay_names=NULL,
 se_names=c("A", "B"),
 startN=2,
 verbose=TRUE,
 ...)
{
   # ensure se_names has at least 2 unique values
   if (length(unique(se_names)) < 2) {
      se_names <- c("A", "B")
   }
   # retrieve existing values upfront
   rowData_colnames_SE1 <- colnames(SummarizedExperiment::rowData(SE1))
   rowData_colnames_SE2 <- colnames(SummarizedExperiment::rowData(SE2))
   colData_colnames_SE1 <- colnames(SummarizedExperiment::colData(SE1));
   colData_colnames_SE2 <- colnames(SummarizedExperiment::colData(SE2));

   # validate rowname1,rowname2
   #
   # First fill rownames where relevant
   if ("rownames" %in% rowname1) {
      if (!"rownames" %in% rowData_colnames_SE1) {
         SummarizedExperiment::rowData(SE1)$rownames <- rownames(SE1);
         rowData_colnames_SE1 <- c(rowData_colnames_SE1, "rownames");
      }
   }
   if ("rownames" %in% rowname2) {
      if (!"rownames" %in% rowData_colnames_SE2) {
         SummarizedExperiment::rowData(SE2)$rownames <- rownames(SE2);
         rowData_colnames_SE2 <- c(rowData_colnames_SE2, "rownames");
      }
   }
   #
   # when rowname1 or rowname2 contain multiple values
   # concatenate into "merge_key"
   if (length(rowname1) > 1 || length(rowname2) > 1) {
      SummarizedExperiment::rowData(SE1)$merge_key <- jamba::pasteByRow(
         data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            SummarizedExperiment::rowData(SE1)[,rowname1, drop=FALSE]));
      SummarizedExperiment::rowData(SE2)$merge_key <- jamba::pasteByRow(
         data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            SummarizedExperiment::rowData(SE2)[,rowname2, drop=FALSE]));
      rowname1 <- "merge_key";
      rowname2 <- "merge_key";
   }
   # when rowname1 and rowname2 differ, create a column
   if (!all(rowname1 == rowname2)) {
      SummarizedExperiment::rowData(SE1)$merge_key <- (
         SummarizedExperiment::rowData(SE1)[[rowname1]]);
      SummarizedExperiment::rowData(SE2)$merge_key <- (
         SummarizedExperiment::rowData(SE2)[[rowname2]]);
   }

   # Prepare to merge these two rowData
   if (TRUE %in% rowData_colnames_intersect) {
      rowData_colnames <- intersect(
         rowData_colnames_SE1,
         rowData_colnames_SE2);
      # verbose output
      if (verbose) {
         if (length(rowData_colnames_SE1) < length(rowData_colnames)) {
            jamba::printDebug("merge_proteomics_se(): ",
               "rowData colnames dropped from SE1:",
               setdiff(rowData_colnames_SE1, rowData_colnames))
         }
         if (length(rowData_colnames_SE2) < length(rowData_colnames)) {
            jamba::printDebug("merge_proteomics_se(): ",
               "rowData colnames dropped from SE2:",
               setdiff(rowData_colnames_SE2, rowData_colnames))
         }
         jamba::printDebug("merge_proteomics_se(): ",
            "rowData colnames retained:",
            rowData_colnames)
      }
      # trim rowData colnames
      SummarizedExperiment::rowData(SE1) <- (
         SummarizedExperiment::rowData(SE1)[,rowData_colnames])
      SummarizedExperiment::rowData(SE2) <- (
         SummarizedExperiment::rowData(SE2)[,rowData_colnames])
   } else {
      stop("rowData_colnames_intersect=FALSE is not yet implemented.")
   }

   # Prepare to merge these two colData
   if (TRUE %in% colData_colnames_intersect) {
      colData_colnames <- intersect(
         colData_colnames_SE1,
         colData_colnames_SE2)
      # verbose output
      if (verbose) {
         if (length(colData_colnames_SE1) < length(colData_colnames)) {
            jamba::printDebug("merge_proteomics_se(): ",
               "colData colnames dropped from SE1:",
               setdiff(colData_colnames_SE1, colData_colnames))
         }
         if (length(colData_colnames_SE2) < length(colData_colnames)) {
            jamba::printDebug("merge_proteomics_se(): ",
               "colData colnames dropped from SE2:",
               setdiff(colData_colnames_SE2, colData_colnames))
         }
         jamba::printDebug("merge_proteomics_se(): ",
            "colData colnames retained:",
            colData_colnames)
      }
      # trim rowData colnames
      SummarizedExperiment::colData(SE1) <- (
         SummarizedExperiment::colData(SE1)[, colData_colnames])
      SummarizedExperiment::colData(SE2) <- (
         SummarizedExperiment::colData(SE2)[, colData_colnames])
   } else {
      stop("colData_colnames_intersect=FALSE is not yet implemented.")
   }

   # Replace empty rowname1 value with rownames(SE1)
   rowData(SE1)[[rowname1]] <- ifelse(
      nchar(SummarizedExperiment::rowData(SE1)[[rowname1]]) > 0,
      SummarizedExperiment::rowData(SE1)[[rowname1]],
      rownames(SE1))
   rownames(SE1) <- jamba::makeNames(
      SummarizedExperiment::rowData(SE1)[[rowname1]],
      renameFirst=FALSE,
      startN=startN,
      ...)
   # Replace empty rowname2 value with rownames(SE2)
   rowData(SE2)[[rowname2]] <- ifelse(
      nchar(SummarizedExperiment::rowData(SE2)[[rowname2]]) > 0,
      SummarizedExperiment::rowData(SE2)[[rowname2]],
      rownames(SE2))
   rownames(SE2) <- jamba::makeNames(
      SummarizedExperiment::rowData(SE2)[[rowname2]],
      renameFirst=FALSE,
      startN=startN,
      ...)

   # prepare to merge unified rownames
   rownames_union <- jamba::mixedSort(union(
      rownames(SE1),
      rownames(SE2)))

   # expand rowData_colnames_unique
   if (any(rowData_colnames_unique %in% rowData_colnames)) {
      rowData_colnames_unique_to <- paste0(rowData_colnames_unique,
         "_se",
         se_names[2]);
      # rename columns in SE2
      SummarizedExperiment::rowData(SE2) <- jamba::renameColumn(
         x=SummarizedExperiment::rowData(SE2),
         from=rowData_colnames_unique,
         to=rowData_colnames_unique_to);
      rowData_colnames_use2 <- colnames(SummarizedExperiment::rowData(SE2))
      if (verbose && !all(rowData_colnames_use2 %in% rowData_colnames)) {
         jamba::printDebug("merge_proteomics_se(): ",
            "Renamed columns in rowData(SE2) to keep unique values:",
            setdiff(rowData_colnames_use2, rowData_colnames))
      }
   } else {
      rowData_colnames_use2 <- rowData_colnames
   }

   # merge rowData by creating empty shell
   # then fill with relevant rows and columns from rowData(SE1), rowData(SE2)
   rowData12_row <- head(
      SummarizedExperiment::rowData(SE1), 1)
   if (!identical(rowData_colnames, rowData_colnames_use2)) {
      # append any renamed columns
      rowData12_row[, rowData_colnames_use2] <- head(
         SummarizedExperiment::rowData(SE2), 1)
      # re-order using original column order so versioned columns are
      # adjacent
      rowData_colnames_order <- jamba::provigrep(
         paste0("^", rowData_colnames, "(|_se", se_names[2], ")$"),
         colnames(rowData12_row))
      rowData12_row <- rowData12_row[, rowData_colnames_order, drop=FALSE]
      if (verbose) {
         jamba::printDebug("merge_proteomics_se(): ",
            "Re-ordered columns to retain the original order of annotations:",
            rowData_colnames_order)
      }
   }
   # convert all values to NA and retain the column class
   for (icol in colnames(rowData12_row)) {
      rowData12_row[1, icol] <- NA;
   }
   # expand to appropriate number of rows
   rowData12 <- rowData12_row[rep(1, length(rownames_union)),, drop=FALSE]
   # assign fixed rownames
   rownames(rowData12) <- rownames_union

   # assign rowData values into this empty shell
   # assign SE2 first, so that SE1 values take priority for any overlaps
   rowData12[rownames(SE2), rowData_colnames_use2] <- (
      SummarizedExperiment::rowData(SE2)[, rowData_colnames_use2, drop=FALSE])
   rowData12[rownames(SE1), rowData_colnames] <- (
      SummarizedExperiment::rowData(SE1)[, rowData_colnames, drop=FALSE])
   # define incidence matrix to show which rows occurred in which Set
   rowData12[, se_names[1]] <- (rownames_union %in% rownames(SE1)) * 1
   rowData12[, se_names[2]] <- (rownames_union %in% rownames(SE2)) * 1


   # confirm that colnames(SE1) and colnames(SE2) are distinct
   while (any(colnames(SE1) %in% colnames(SE2))) {
      colnames_SE1 <- colnames(SE1);
      colnames_SE2 <- colnames(SE2);
      colnames_SE12 <- jamba::makeNames(
         c(colnames_SE1, colnames_SE2),
         renameFirst=FALSE,
         startN=startN,
         ...);
      colnames(SE1) <- head(colnames_SE12,
         ncol(SE1));
      colnames(SE2) <- tail(colnames_SE12,
         ncol(SE2));
      if (verbose) {
         if (length(colData_colnames_SE1) < length(colData_colnames)) {
            jamba::printDebug("merge_proteomics_se(): ",
               "colnames in SE1 SE2 were renamed to be unique.")
            jamba::printDebug("merge_proteomics_se(): ",
               "colnames(SE1):", colnames(SE1));
            jamba::printDebug("merge_proteomics_se(): ",
               "colnames(SE2):", colnames(SE2));
         }
      }
   }
   # unified colnames
   colnames_union <- union(
      colnames(SE1),
      colnames(SE2))

   # merge colData12 by rbind()
   colData12 <- rbind(
      SummarizedExperiment::colData(SE1),
      SummarizedExperiment::colData(SE2))[colnames_union, , drop=FALSE]
   # colData12$Set <- rep(se_names[1:2], c(ncol(SE1), ncol(SE2)));
   colData12[, se_names[1]] <- (colnames_union %in% colnames(SE1)) * 1
   colData12[, se_names[2]] <- (colnames_union %in% colnames(SE2)) * 1

   # prepare merged assays
   if (length(assay_names) == 0) {
      assay_names <- assayNames(SE1)
   }
   assay_names <- jamba::nameVector(assay_names);
   row_match1 <- match(rownames_union, rownames(SE1));
   row_match2 <- match(rownames_union, rownames(SE2));
   assay_list <- lapply(assay_names, function(iassay){
      assay1 <- assays(SE1[row_match1, , drop=FALSE])[[iassay]];
      assay2 <- assays(SE2[row_match2, , drop=FALSE])[[iassay]];
      m <- cbind(
         assays(SE1[row_match1, , drop=FALSE])[[iassay]],
         assays(SE2[row_match2, , drop=FALSE])[[iassay]])
      rownames(m) <- rownames_union;
      colnames(m) <- colnames_union;
      m;
   })

   # create combined SummarizedExperiment
   SE12 <- SummarizedExperiment(
      assays=assay_list,
      colData=colData12,
      rowData=rowData12)
   return(SE12)
}
