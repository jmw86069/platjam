
#' Import metabolomics data from NIEHS file formats
#'
#' Import metabolomics data from NIEHS file formats
#'
#' This import function is specific to NIEHS file formats produced
#' from their defined analysis workflow. The files typically include
#' `"df_pos"` for positive ionization, and `"df_neg"` for negative
#' ionization:
#'
#' ## Sample Metadata
#'
#' * `"[project_code]_NIEHS_MCF_metadata.txt"`: tab-delimited text file
#' which contains sample annotations.
#'
#' ## Positive Ionization Files
#'
#' * `"df_pos.datamatrix.cleaned.txt"`: Tab-delimited text file
#' containing peak areas.
#' Features are processed and cleaned by MCF for quality.
#' * `"df_pos.datamatrix.cleaned.log10.txt"`: As above, but log10 transformed.
#' * `"df_pos.datamatrix.cleaned.rowsum.txt"`: As above but using
#' row sum peak areas.
#' * `"df_pos.annotation.cleaned.txt"`: annotation of each measured metabolite.
#'
#' ## Negative Ionization Files
#'
#' * `"df_neg.datamatrix.cleaned.txt"`: Tab-delimited text file
#' containing peak areas.
#' Features are processed and cleaned by MCF for quality.
#' * `"df_neg.datamatrix.cleaned.log10.txt"`: As above, but log10 transformed.
#' * `"df_neg.datamatrix.cleaned.rowsum.txt"`: As above but using
#' row sum peak areas.
#' * `"df_neg.annotation.cleaned.txt"`: annotation of each measured metabolite.
#'
#' @family jam import functions
#'
#' @return `list` of `SummarizedExperiment` objects, where the `list`
#'    is defined by the type of ionization ("df_pos", "df_neg"), and
#'    the type of data ("cleaned") in the data filenames.
#'    Typically the result includes:
#'    * `"df_pos_cleaned"`
#'    * `"df_neg_cleaned"`
#'
#'    If only one ionization is provided, only one entry will be returned.
#'
#'    For each `SummarizedExperiment` object:
#'    * `rowData` represents metabolite annotations
#'    * `colData` represents sample annotations, optionally including
#'    annotations via a `data.frame` supplied as `curation_txt`.
#'    * the slot `"metadata"` is a `list` with the following:
#'
#'       * `isample_use`: the subset of `colnames(se)` for which there
#'       was sample metadata found in the metadata file. Some control
#'       samples may not match the full metadata, and will be ignored
#'       when using `isamples_Use`.
#'       * `irows_use`: all `rownames(se)` for all measured metabolites.
#'       * `irows_clean`: the `rownames(se)` for measurement with no
#'       annotation in the column `"flag_guidance"`.
#'       * `irows_flagged`: the `rownames(se)` for measurements with
#'       some non-empty annotation in the column `"flag_guidance"`.
#'
#' @param data_path `character` path which contains files as named in
#'    Details.
#' @param shared_samples_only `logical` indicating whether to retain
#'    only those samples which are shared across positive and negative
#'    ionization data files. This step also typically removed extraneous
#'    technical QC samples which may not have been included in both
#'    the positive and negative ionization data.
#' @param filter_sample_type `character` with optional subset of
#'    sample types to retain. The most useful is
#'    `filter_sample_types="sample"` which will retain only biological
#'    samples, and will drop all other samples. The commonly observed
#'    values in the column `"CORE_SampleType"`:
#'    * `"sample"`: biological sample
#'    * `"QC_curve"`: calibration curve sample
#'    * `"blank_system"`: negative control blank sample
#'    * `"<NA>"`: empty values where data provided a sample which was
#'    not described in the associated project metadata file. This
#'    typically occurs only for system quality control checks,
#'    such as "SystemSuitability", "wash", "AqX_blank", and "AqX_sample".
#' @param curation_txt `data.frame` passed to `curate_se_colData()`
#'    in order to include sample annotations. The default uses
#'    identifiers from `colnames(se)` for each `SummarizedExperiment`
#'    object.
#'    column headers found in the annotation metadata file.
#'    If `curation_txt` is not supplied, then values will be split into
#'    columns by `_` underscore or `" "` whitespace characters.
#' @param drop_na_columns `logivcal` indicating whether to drop columns
#'    in `colData` or `rowData` when all values are `NA`, `"NA"`,
#'    `"not applicable"`, or `"none"`, defined in `drop_na_values`.
#'
#'    Note that columns with only one non-na value in all fields will
#'    be removed from the `colData(se)` and stored as metadata as
#'    a single `character` vector accessible via:
#'    `se@metadata$colData_values`
#' @param drop_na_values `character` vector of values considered to
#'    be "na" when `drop_na_columns=TRUE`.
#' @param simplify_singlet_columns `logical` indicating whether to remove
#'    `colData(se)` columns with only one value, instead storing the
#'    name and value in metadata as a character vector.
#'    This option is enabled by default, and simplifies the resulting
#'    `colData(se)` so that it only includes columns with two or
#'    more unique values.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' data_path <- path.expand("~/Projects/Rider/metabolomics_jul2023/data");
#' se_list <- import_metabolomics_niehs(data_path);
#'
#' @export
import_metabolomics_niehs <- function
(data_path,
 shared_samples_only=TRUE,
 filter_sample_type=NULL,
 curation_txt=NULL,
 drop_na_columns=TRUE,
 drop_na_values=c(NA,
    "NA",
    "not applicable",
    "not specified",
    "none"),
 simplify_singlet_columns=TRUE,
 verbose=FALSE,
 ...)
{
   # validate data_path
   if (length(data_path) == 0 || nchar(data_path) == 0) {
      stop("No data_path was supplied.")
   }
   if (!dir.exists(data_path)) {
      stop(paste0("The data_path was not accessible: '", data_path, "'."))
   }

   # validate curation_txt
   if (length(curation_txt) > 0) {
      if (length(curation_txt) == 1 &&
            is.atomic(curation_txt)) {
         if (file.exists(curation_txt)) {
            if (verbose) {
               jamba::printDebug("import_metabolomics_niehs(): ",
                  "Importing curation_txt from file '", curation_txt, "'.")
            }
            curation_txt <- data.table::fread(curation_txt,
               data.table=FALSE)
         } else {
            stop("curation_txt must be either data.frame or filename.")
         }
      }
   }

   # find metadata file
   metadata_file <- list.files(path=data_path,
      full.names=TRUE,
      ignore.case=TRUE,
      pattern="^[a-zA-Z0-9].*_NIEHS_MCF_Metadata.txt$")
   if (length(metadata_file) == 0) {
      stop(paste0(
         "No metadata file was found using suffix '",
         "_NIEHS_MCF_Metadata.txt", "'."))
   }
   if (length(metadata_file) > 1) {
      stop(paste0(
         "There should be only one metadata file during import, but ",
         length(metadata_file), " were found."))
   }
   # load sample metadata
   metadata_df <- data.table::fread(metadata_file,
      data.table=FALSE)
   # 0.0.67.900 - fix issue where CORE_Filename may be present multiple times
   metadata_df <- subset(metadata_df, !duplicated(CORE_Filename))
   # now rownames should be unique
   rownames(metadata_df) <- metadata_df$CORE_Filename;

   # find data files
   data_files <- list.files(path=data_path,
      full.names=TRUE,
      ignore.case=TRUE,
      pattern="^df_([a-z]+)[.]([a-z]+)[.](.*)[.]txt$")
   data_files_df <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      jamba::rbindList(strsplit(
            gsub("[.]txt$", "", basename(data_files)),
         "[.]")))
   rownames(data_files_df) <- basename(data_files);
   colnames(data_files_df) <- c("Ionization",
      "DataType",
      "RowType",
      "Signal")

   data_files_list <- split(data_files_df,
      jamba::pasteByRowOrdered(
         data_files_df[,c("Ionization", "RowType")]))
   # iterate each subset of files
   se_list <- lapply(jamba::nameVectorN(data_files_list), function(ilist){
      if (verbose) {
         jamba::printDebug("import_metabolomics_niehs(): ",
            "Importing file subset '", ilist, "'.")
      }
      idf <- data_files_list[[ilist]];
      # import annotation
      if ("annotation" %in% idf$DataType) {
         idf_sub <- head(subset(idf, DataType %in% "annotation"), 1)
         anno_df <- data.table::fread(
            file.path(data_path, rownames(idf_sub)),
            data.table=FALSE)
         fi_colname <- head(jamba::vigrep("^feature.index$",
            colnames(anno_df)), 1)
         if (length(fi_colname) == 1) {
            rownames(anno_df) <- anno_df[[fi_colname]];
         } else {
            jamba::printDebug("import_metabolomics_niehs(): ",
               "No annotation column 'feature_index' exists for subset '",
               ilist,
               "'. Skipping this subset.",
               fgText=c("darkorange", "red"))
            return(NULL)
         }
      } else {
         jamba::printDebug("import_metabolomics_niehs(): ",
            "No annotation file exists for subset '",
            ilist,
            "'. Skipping this subset.",
            fgText=c("darkorange", "red"))
         return(NULL)
      }
      if ("datamatrix" %in% idf$DataType) {
         idf_sub <- subset(idf, DataType %in% "datamatrix")
         idf_sub$Signal <- gsub("^$", "data", idf_sub$Signal)
         idf_sub <- jamba::mixedSortDF(idf_sub, byCols="Signal")
         idf_files <- jamba::nameVector(rownames(idf_sub), idf_sub$Signal);
         assay_list <- lapply(idf_files, function(idf_file) {
            jdf <- data.table::fread(
               file.path(data_path, idf_file),
               data.table=FALSE)
            sample_colname <- head(
               jamba::vigrep("^SampleID$", colnames(jdf)),
               1)
            if (length(sample_colname) == 1) {
               rownames(jdf) <- jdf[[sample_colname]];
               use_colnames <- setdiff(colnames(jdf), sample_colname)
               jmatrix <- t(jdf[, use_colnames, drop=FALSE])
               return(jmatrix)
            }
            jamba::printDebug("import_metabolomics_niehs(): ",
               "No measurement column 'SampleID' exists for subset '",
               ilist,
               "' and file '",
               ifile,
               "'. Skipping this file.",
               fgText=c("darkorange", "red"))
            return(NULL)
         })
         assay_rows <- rownames(assay_list[[1]]);
         assay_columns <- colnames(assay_list[[1]]);
      } else {
         jamba::printDebug("import_metabolomics_niehs(): ",
            "No datamatrix files exists for subset '",
            ilist,
            "'. Skipping this subset.",
            fgText=c("darkorange", "red"))
         return(NULL)
      }

      # colData: align annotation with assay data
      column_match <- match(assay_columns,
         metadata_df$CORE_Filename)
      colData_df <- metadata_df[column_match, , drop=FALSE];
      rownames(colData_df) <- assay_columns;
      # add common SampleID without _p_ (positive) or _n_ (negative)
      colData_df$SampleID <- gsub("_[np]_", "_", assay_columns);

      # rowData: align annotation with assay data
      row_match <- match(assay_rows,
         anno_df[[fi_colname]]);
      rowData_df <- anno_df[row_match, , drop=FALSE];
      rownames(rowData_df) <- assay_rows;

      # define usable sample colnames
      isamples_use <- rownames(
         subset(colData_df, CORE_SampleType %in% "sample"))

      # define usable rownames
      irows_use <- assay_rows;
      # note it is unclear whether NA is different than ""
      # however NA seems to indicate metabolites with no search results,
      # while "" seems to indicate metabolites with some database result.
      irows_clean <- rownames(
         subset(rowData_df, flag_guidance %in% c("", NA)))
      irows_flagged <- setdiff(irows_use, irows_clean);

      # note that some combinations "Formula", "Calc. MW" are duplicated
      # (skip for now)
      if (FALSE) {
         formula_calcmw <- jamba::pasteByRow(
            rowData_df[,c("Formula", "Calc. MW"), drop=FALSE])
         formula_calcmw_dupe <- jamba::tcount(formula_calcmw, 2)
         head(formula_calcmw_dupe, 10)
         rowData_dupe_df <- subset(rowData_df,
            formula_calcmw %in% names(formula_calcmw_dupe)[2])
         rowData_dupe_df
      }

      # optionally remove columns whose values are all "NA"/"not applicable"
      # colData
      # colData_colnames_uvalues <- lapply(jamba::nameVectorN(colData_df), function(icol){
      #    vals <- colData_df[[icol]];
      #    if (TRUE %in% drop_na_columns && length(drop_na_values) > 0) {
      #       vals[vals %in% drop_na_values] <- NA;
      #    }
      #    uvals <- jamba::rmNA(unique(vals));
      #    uvals
      # })
      # if (TRUE %in% drop_na_columns && length(drop_na_values) > 0) {
      #    colData_colnames_nvalues <- lengths(colData_colnames_uvalues)
      #    colData_metadata <- NULL;
      #    colData_singlets <- (colData_colnames_nvalues == 1);
      #    if (any(colData_singlets)) {
      #       colData_metadata <- unlist(
      #          colData_colnames_uvalues[colData_singlets])
      #    }
      #    colData_colnames_keep <- colnames(colData_df)[
      #       colData_colnames_nvalues > 1]
      #    # colData_df <- colData_df[, colData_colnames_keep, drop=FALSE];
      #
      #    # rowData is skipped for now.
      #    # All values are unique except "mode" which indicates
      #    # "positive" or "negative" and should probably be kept as-is.
      # } else {
      #    colData_colnames_keep <- colnames(colData_df)
      # }

      # assemble SummarizedExperiment
      se <- SummarizedExperiment::SummarizedExperiment(
         assays=assay_list,
         rowData=rowData_df,
         colData=colData_df,
         metadata=list(
            isamples_use=isamples_use,
            irows_use=irows_use,
            irows_clean=irows_clean,
            irows_flagged=irows_flagged
            # colData_values=colData_colnames_uvalues,
            # colData_colnames_keep=colData_colnames_keep
         ))

      # optionally apply curation_txt
      if (length(curation_txt) > 0) {
         se <- curate_se_colData(se,
            df=curation_txt,
            verbose=verbose,
            indent=3,
            ...)
         # jamba::printDebug("colnames(colData(se)):",
         #    colnames(SummarizedExperiment::colData(se)))
      }

      return(se)
   })

   # revise colData for each entry
   # to remove columns with no values,
   # and simplify columns with only one unique value.
   se_list <- lapply(jamba::nameVectorN(se_list), function(se_name) {
      se <- se_list[[se_name]];
      colData_colnames <- colnames(SummarizedExperiment::colData(se));

      # define the unique non-NA values in each column
      colData_colnames_uvalues <- lapply(jamba::nameVector(colData_colnames),
         function(icol){
         vals <- SummarizedExperiment::colData(se)[[icol]];
         if (TRUE %in% drop_na_columns && length(drop_na_values) > 0) {
            vals[vals %in% drop_na_values] <- NA;
         }
         if (all(is.na(vals))) {
            uvals <- NULL
         } else {
            uvals <- jamba::rmNA(unique(vals));
         }
         uvals
      })
      # optionally remove columns whose values are all "NA"/"not applicable"
      colData_metadata <- NULL;
      colData_colnames_nvalues <- lengths(colData_colnames_uvalues)
      if (any(colData_colnames_nvalues < 2)) {
         if (TRUE %in% simplify_singlet_columns) {
            colData_singlets <- (colData_colnames_nvalues == 1);
            if (any(colData_singlets)) {
               colData_metadata <- unlist(
                  colData_colnames_uvalues[colData_singlets])
            }
            colData_colnames_keep <- colData_colnames[
               colData_colnames_nvalues > 1]
         } else {
            colData_colnames_keep <- colData_colnames[
               colData_colnames_nvalues > 0]
         }
      } else {
         colData_colnames_keep <- colData_colnames
      }
      se@metadata$colData_values <- colData_colnames_uvalues;
      se@metadata$colData_colnames_keep <- colData_colnames_keep;
      se;
   })

   # apply colData_colnames_keep across all se_list entries for consistency
   colData_colnames_list <- lapply(se_list, function(se){
      se@metadata$colData_colnames_keep
   })
   colData_colnames <- Reduce("union", colData_colnames_list);

   # keep columns with multiple values in any se_list entry
   # move columns with one non-na entry into metadata
   se_list <- lapply(jamba::nameVectorN(se_list), function(se_name){
      se <- se_list[[se_name]];
      # remove this temporary metadata value
      se@metadata$colData_colnames_keep <- NULL;
      # define colnames to keep
      colData_colnames_use <- intersect(
         colnames(SummarizedExperiment::colData(se)),
         colData_colnames)
      # define colnames to drop
      colData_colnames_drop <- setdiff(
         colnames(SummarizedExperiment::colData(se)),
         colData_colnames)
      if (length(colData_colnames_drop) == 0) {
         se@metadata$colData_values <- NULL;
      } else {
         se@metadata$colData_values <- jamba::rmNULL(
            se@metadata$colData_values[colData_colnames_drop])
         if (all(lengths(se@metadata$colData_values) == 1)) {
            se@metadata$colData_values <- unlist(se@metadata$colData_values);
         }
         if (verbose) {
            jamba::printDebug("import_metabolomics_niehs(): ",
               "Reducing colData in '",
               se_name,
               "' from ",
               ncol(SummarizedExperiment::colData(se)),
               " columns to ",
               length(colData_colnames_use))
         }
         SummarizedExperiment::colData(se) <- SummarizedExperiment::colData(se)[,
            colData_colnames_use, drop=FALSE]
      }
      return(se);
   })

   return(se_list)
}
