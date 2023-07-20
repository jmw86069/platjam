
#' Curate Summarized Experiment colData
#'
#' Apply curation to colData in a SummarizedExperiment object
#'
#' Given a SummarizedExperiment object, this function is intended
#' to augment the `SummarizedExperiment::colData()` annotation associated
#' with columns, which are typically biological or experimental
#' samples.
#' Measurements within each sample are typically stored as rows.
#'
#' A convenient wrapper to `curate_to_df_by_pattern()`, which applies
#' the result directly to `SummarizedExperiment::colData()` which is
#' stored as a `S4Vectors::DataFrame-class`.
#'
#' Note that colnames present in both `colData(se)` and `df` will
#' take the value from `df` as replacement, including the presence of `NA`
#' values.
#'
#' ## About pattern matching
#'
#' The patterns are used to match identifiers using regular expressions,
#' and the argument `warn_multimatch=TRUE` (default) will print a
#' warning when one pattern matches two or more identifiers.
#' It may be intended, or may indicate that some patterns are not
#' specific enough to match only one intended identifier.
#'
#' For example `pattern="sample_3"` will match identifiers:
#' `c("one_sample_3", "two_sample_3", "one_sample_31")`.
#'
#' To overcome this type of issue, use regular expressions to
#' limit matching to the end, for example `pattern="sample_3$"`
#' will only match `c("one_sample_3", "two_sample_3")` and
#' will not match `"one_sample_31"`.
#'
#' It can be helpful to name the pattern column `"Pattern"` so that
#' the pattern used is clearly defined in the output
#' `colData(se)`, and can be compared to the intended identifiers.
#'
#' @return `SummarizedExperiment::SummarizedExperiment` object.
#'    * When `subset_se=FALSE` (default), the output will contain
#'    the same dimensions and column order as the input `se`.
#'    * When `subset_se=TRUE` the output object may contain fewer columns
#'    based upon the number of identifiers that matched the patterns
#'    supplied in `df`.
#'
#' @family jam utility functions
#'
#' @param se `SummarizedExperiment` object.
#' @param df `data.frame` (or equivalent) which contains columns of
#'    data annotation to be applied.
#'    The first column is assumed to be the column used for
#'    patterns to be matched with identifiers in the `se` object.
#'    The pattern column can be defined with `pattern_colname`.
#' @param pattern_colname `character` value indicating which
#'    column in `df` contains patterns to be matched with identifiers
#'    in `se`. The default uses the first column in `df`.
#'    This value is passed to `curate_to_df_by_pattern()`.
#' @param group_colname `character` or `NULL` (default) indicating
#'    which column(s) represent experimental groups, used
#'    only to create a corresponding column with unique label
#'    for each entry. When `NULL` no action is taken, which is default.
#' @param id_colname `character` used only when `group_colname` is
#'    defined and present in `colnames(df)`, used to create a
#'    unique label for each row in `colData(se)`.
#'    By default `group_colname=NULL` so no action is taken.
#' @param use `character` string indicating the data to use as
#'    the identifiers when applying curation logic.
#'    The default is to use `colnames(se)`, however it can use one
#'    or more columns from `SummarizedExperiment::colData(se)`.
#'    Some options are described below:
#'    * `"colnames"`: uses `colnames(se)`, which should be equivalent
#'    to using `rownames(SummarizedExperiment::colData(se))`.
#'    * `"rownames"`: uses `rownames(SummarizedExperiment::colData(se))`,
#'    which as stated above should be equivalent to `colnames(se)`.
#'    * one or more `character` values that match `colnames(colData(se))`.
#' @param use_delim `character` string used as a delimiter when
#'    `use` is supplied as a vector with multiple colnames.
#'    The values in each column are concatenated using this delimiter,
#'    by calling `jamba::pasteByRow()`.
#' @param subset_se `logical` indicating whether the `se` object columns
#'    be subset when not all identifiers matched the patterns in `df`.
#'    * When `subset_se=FALSE` any entries in `se` for which
#'    the identifier did not match the pattern in `df`,
#'    the corresponding rows of `SummarizedExperiment::colData()`
#'    will contain `NA` values.
#'    * When `subset=TRUE` any entries in `se` for which
#'    the identifier did not match the pattern in `df` will
#'    be removed from the `se` object. This option is sometimes
#'    a convenient way to subset a large data to use only
#'    user-defined samples.
#' @param warn_multimatch `logical` indicating whether to print a warning
#'    when any one pattern matches two or more identifiers.
#'    Sometimes this behavior is intended, however it may indicate
#'    that the patterns are not specific enough to match one unique
#'    identifier. See Details.
#' @param indent `numeric` value used when `verbose=TRUE`, passed to
#'    `jamba::printDebug()`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `curate_to_df_by_pattern()`.
#'
#' @export
curate_se_colData <- function
(se,
 df,
 pattern_colname=head(colnames(df), 1),
 group_colname=NULL,
 id_colname="Label",
 use=c("colnames"),
 use_delim="_",
 subset_se=FALSE,
 warn_multimatch=TRUE,
 indent=0,
 verbose=TRUE,
 ...)
{
   # handle argument 'use'
   if (length(use) == 0) {
      stop("The argument 'use' must be supplied, but was empty.")
   }
   # convert rownames to colnames, take only unique values
   use <- unique(gsub("^rownames$", "colnames", use));
   if ("colnames" %in% use) {
      if (length(use) > 1) {
         stop(paste("Argument 'use' must only define 'colnames'",
            "or colnames present in colnames(colData(se)), not both."))
      }
      identifiers <- colnames(se);
      input_colname <- jamba::ucfirst("colnames");
   } else {
      colData_colnames <- colnames(SummarizedExperiment::colData(se));
      if (!all(use %in% colData_colnames)) {
         missing_names <- setdiff(use, colData_colnames);
         stop(paste0(
            "All values in 'use' must be present in colnames(colData(se)), ",
            "and some were not: ",
            jamba::cPaste(paste0("'", missing_names, "'"))
            ))
      }
      identifiers <- jamba::pasteByRow(data.frame(stringsAsFactors=FALSE,
         SummarizedExperiment::colData(se)[, use, drop=FALSE]),
         sep=use_delim)
      input_colname <- jamba::ucfirst(jamba::cPaste(use, sep=use_delim));
   }

   # call curation function
   new_colData <- curate_to_df_by_pattern(
      x=identifiers,
      pattern_colname=pattern_colname,
      input_colname=input_colname,
      group_colname=group_colname,
      id_colname=id_colname,
      df=df,
      verbose=FALSE,
      ...,
      order_priority="x")

   # check for duplicate values in pattern_colname, which may indicate
   # the pattern was not specific enough to match only one entry.
   if (nrow(new_colData) > 0) {
      pattern_tc <- jamba::tcount(new_colData[[pattern_colname]], 2)
      if (length(pattern_tc) > 0) {
         if (verbose) {
            jamba::printDebug("curate_se_colData(): ",
               "Warning: ",
               "patterns matched multiple identifiers: '",
               head(names(pattern_tc), 1), "' was matched ",
               head(pattern_tc, 1), " times.",
               indent=indent,
               fgText=c("darkorange", "orangered"))
            jamba::printDebug("curate_se_colData(): ",
               "Warning: ",
               "Verify the patterns match the intended, specific identifiers.",
               indent=indent,
               fgText=c("darkorange", "orangered"))
         }
      }
   } else {
      if (verbose) {
         jamba::printDebug("curate_se_colData(): ",
            "Warning: ",
            "No patterns were matched.",
            indent=indent,
            fgText=c("darkorange", "orangered"))
      }
   }
   # if nrow() == 0 there is no match


   # handle missing data
   if (nrow(new_colData) < ncol(se)) {
      num_missing <- ncol(se) - nrow(new_colData);
      if (FALSE %in% subset_se) {
         if (verbose) {
            jamba::printDebug("curate_se_colData(): ",
               indent=indent,
               "Introducing NA values for ",
               jamba::formatInt(num_missing),
               " missing entries in se.")
         }
         id_match <- match(identifiers,
            new_colData[[input_colname]])
         new_colData <- new_colData[id_match, , drop=FALSE]
         rownames(new_colData) <- identifiers
      } else {
         if (verbose) {
            jamba::printDebug("curate_se_colData(): ",
               indent=indent,
               "Subsetting se to remove ",
               jamba::formatInt(num_missing),
               " missing entries.")
         }
         colData_match <- match(new_colData[[input_colname]],
            identifiers);
         se <- se[, colData_match];
         # subset se
      }
   } else {
      if (verbose) {
         jamba::printDebug("curate_se_colData(): ",
            indent=indent,
            "All ",
            jamba::formatInt(ncol(se)),
            " identifiers were matched.")
      }
   }

   # apply new_colData into colData(se)
   new_colnames <- setdiff(colnames(new_colData),
      "Colnames")
   if (verbose) {
      jamba::printDebug("curate_se_colData(): ",
         indent=indent,
         "Updating colData colnames: ",
         paste0("'", new_colnames, "'"))
   }
   SummarizedExperiment::colData(se)[, new_colnames] <- (
      new_colData[, new_colnames, drop=FALSE]);

   # return the modified se object
   return(se);
}
