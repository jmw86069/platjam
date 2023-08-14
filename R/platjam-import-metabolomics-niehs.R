
#' Import metabolomics data from NIEHS file formats
#'
#' Import metabolomics data from NIEHS file formats
#'
#' This import function is specific to NIEHS file formats produced
#' from their defined analysis workflow. The files typically include
#' `"df_pos"` for positive ionization, and `"df_neg"` for negative
#' ionization.
#'
#' Optionally, when the full data processing file is present,
#' it will be imported alongside the cleaned data described above.
#' The full data processing imports detailed compound measurement
#' data, and is expected in one of two formats in the `data_path` folder:
#' 1. Files `"compounds_pos.txt"` and/or `"compounds_neg.txt"`, or
#' 2. `"1_DataProcessed.zip"` which is expected to contain
#' files `"compounds_pos.txt"` and/or `"compounds_neg.txt"` in the archive.
#'
#' The "compounds" data includes important annotations for each measurement,
#' specifically the type of numeric measurement that is supplied
#' by the upstream software. These annotations include whether numeric
#' values were imputed, or measured directly in each sample.
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

   #############################################
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
   # 0.0.68.900 - fix issue where CORE_Filename may be present multiple times
   dupe_filenames <- names(jamba::tcount(metadata_df$CORE_Filename, 2));
   if (length(dupe_filenames) > 0) {
      metadata_df <- subset(metadata_df, !duplicated(CORE_Filename))
      if (verbose) {
         jamba::printDebug("import_metabolomics_niehs(): ",
            c("Warning: ",
               jamba::formatInt(length(dupe_filenames)),
               " CORE_Filename values were duplicated in the metadata.",
               " Retaining non-duplicated rows."),
            sep="",
            fgText=c("darkorange", "firebrick3"))
         print(dupe_filenames);
      }
   }
   # now rownames should be unique
   rownames(metadata_df) <- metadata_df$CORE_Filename;

   #############################################
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
         if (verbose) {
            jamba::printDebug("import_metabolomics_niehs(): ",
               indent=3,
               "Importing file subset '", ilist, "' annotation '",
               rownames(idf_sub), "'")
         }
         fi_colname <- head(jamba::vigrep("^feature.index$",
            colnames(anno_df)), 1)
         if (length(fi_colname) == 1) {
            rownames(anno_df) <- anno_df[[fi_colname]];
         } else {
            jamba::printDebug("import_metabolomics_niehs(): ",
               indent=3,
               "No annotation column 'feature_index' exists for subset '",
               ilist,
               "'. Skipping this subset.",
               fgText=c("darkorange", "red"))
            return(NULL)
         }
      } else {
         jamba::printDebug("import_metabolomics_niehs(): ",
            indent=3,
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
            if (verbose) {
               jamba::printDebug("import_metabolomics_niehs(): ",
                  indent=3,
                  "Importing file subset '", ilist, "'       data '",
                  idf_file, "'")
            }
            jdf <- data.table::fread(
               file.path(data_path, idf_file),
               data.table=FALSE)
            sample_colname <- head(
               jamba::vigrep("^SampleID$", colnames(jdf)),
               1)
            if (length(sample_colname) == 1) {
               # 0.0.68.900 - fix issue where SampleID may be duplicated
               dupe_samples <- names(jamba::tcount(jdf[[sample_colname]], 2));
               if (length(dupe_samples) > 0) {
                  if (verbose) {
                     jamba::printDebug("import_metabolomics_niehs(): ",
                        indent=6,
                        c("Warning: ",
                           jamba::formatInt(length(dupe_samples)),
                           " samples were duplicated in the results.",
                           " Retaining non-duplicated rows."),
                        sep="",
                        fgText=c("darkorange", "firebrick3"))
                     print(dupe_samples);
                  }
                  jdf <- subset(jdf, !duplicated(jdf[[sample_colname]]));
               }
               rownames(jdf) <- jdf[[sample_colname]];
               use_colnames <- setdiff(colnames(jdf), sample_colname)
               jmatrix <- t(jdf[, use_colnames, drop=FALSE])
               return(jmatrix)
            }
            jamba::printDebug("import_metabolomics_niehs(): ",
               indent=3,
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
            indent=3,
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

   #############################################
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

   # 0.0.71.900 - optionally load compound_pos.txt and compound_neg.txt
   dp1zip <- list.files(path=data_path,
      pattern="1_DataProcessed.zip",
      full.names=TRUE)
   compound_files <- list.files(path=data_path,
      pattern="^compound_(pos|neg).txt$",
      full.names=TRUE)
   dp1_compound_pos <- NULL;
   dp1_compound_neg <- NULL;
   #
   if (length(compound_files) < 2 && length(dp1zip) == 1) {
      dp1files <- unzip(dp1zip, list=TRUE)
      dp1_compound_file_pos <- jamba::vigrep("compounds_pos.txt",
         dp1files$Name)
      if (length(dp1_compound_file_pos) > 0) {
         if (!file.exists(dp1_compound_file_pos)) {
            dp1_compound_file_pos <- head(unzip(dp1zip,
               files=dp1_compound_file_pos), 1)
         }
         if (verbose) {
            jamba::printDebug("import_metabolomics_niehs(): ",
               "Loading '", dp1_compound_file_pos,
               "' extracted from ZIP archive.")
         }
         dp1_compound_pos <- data.table::fread(dp1_compound_file_pos,
            data.table=FALSE)
      }
      dp1_compound_file_neg <- jamba::vigrep("compounds_neg.txt",
         dp1files$Name)
      if (length(dp1_compound_file_neg) > 0) {
         if (!file.exists(dp1_compound_file_neg)) {
            dp1_compound_file_neg <- head(unzip(dp1zip,
               files=dp1_compound_file_neg), 1)
         }
         if (verbose) {
            jamba::printDebug("import_metabolomics_niehs(): ",
               "Loading '", dp1_compound_file_neg,
               "' extracted from ZIP archive.")
         }
         dp1_compound_neg <- data.table::fread(dp1_compound_file_neg,
            data.table=FALSE)
      }
   } else if (length(compound_files) > 0) {
      # load pos file
      cpd_file_pos <- jamba::vigrep("_pos", compound_files)
      if (length(cpd_file_pos) == 1) {
         if (verbose) {
            jamba::printDebug("import_metabolomics_niehs(): ",
               "Loading '", cpd_file_pos,
               "'.")
         }
         dp1_compound_pos <- data.table::fread(cpd_file_pos,
            data.table=FALSE)
      }

      # load neg file
      cpd_file_neg <- jamba::vigrep("_neg", compound_files)
      if (length(cpd_file_neg) == 1) {
         if (verbose) {
            jamba::printDebug("import_metabolomics_niehs(): ",
               "Loading '", cpd_file_neg,
               "'.")
         }
         dp1_compound_neg <- data.table::fread(cpd_file_neg,
            data.table=FALSE)
      }
   }
   if (length(dp1_compound_pos) > 0) {
      if (verbose) {
         jamba::printDebug("import_metabolomics_niehs(): ",
            "Merging compounds_pos into se_list.")
      }
      se_list <- process_metab_compounds_file(dp1_compound_pos,
         type="pos",
         se_list=se_list)
   }
   if (length(dp1_compound_neg) > 0) {
      if (verbose) {
         jamba::printDebug("import_metabolomics_niehs(): ",
            "Merging compounds_neg into se_list.")
      }
      se_list <- process_metab_compounds_file(dp1_compound_neg,
         type="neg",
         se_list=se_list)
   }

   # return SE list
   return(se_list)
}

#' Process NIEHS compounds output files (internal)
#'
#' Process NIEHS compounds output files (internal)
#'
#' This function produces `SummarizedExperiment` object with `assays`
#' defined by the `assay_names_extra` patterns below. Each assay name
#' is the lowercase value, with spaces replaced with underscore `"_"`.
#' For example `"Gap Fill Status"` will become assay name `"gap_fill_status"`.
#'
#' @family jam import functions
#'
#' @param x `data.frame` after reading `compound_pos.txt` or `compound_neg.txt`
#' @param nsmall_mz `numeric` argument passed to `format()` to define
#'    the number of decimal positions to retain in the output label
#'    from column `"m/z"`.
#' @param nsmall_rt `numeric` argument passed to `format()` to define
#'    the number of decimal positions to retain in the output label
#'    from column `"RT [min]"`.
#' @param assay_names_extra `character` vector with prefix values expected
#'    in column names, for example where `"Area"` is expected to form a prefix
#'    to column names: `"Area: sample_a", "Area: sample_b"`;
#'    and `"Gap Fill Status"` would form column names
#'    `"Gap Fill Status: sample_a", "Gap Fill Status: sample_b"`.
#' @param mz_colname,rt_colname `character` string to indicate the appropriate
#'    column name containing `"m/z"` (mass/charge) and `"RT [min]"`
#'    (retention time in minutes) as used to format a unique identifier
#'    for each row.
#' @param type `character` string indicating whether data represents
#'    positive `"pos"` or negative `"neg"` ionization. Mainly used
#'    when `se_list` is also provided.
#' @param se_list optional `list` of `SummarizedExperiment` objects,
#'    named with `"pos"` or `"neg"` in the name. When provided, the
#'    `se` output of this function is merged into the `se_list` with
#'    matching `type`, adding `assay_names` into the corresponding
#'    `se_list` object.
#' @param ... additional arguments are ignored.
#'
#' @returns `SummarizedExperiment` object with `assays` defined by matching
#'    column names in `x`. When no column names meet the `assay_names_extra`
#'    criteria, this function returns `NULL`.
#'
#' @export
process_metab_compounds_file <- function
(x,
 nsmall_mz=4,
 nsmall_rt=2,
 assay_names_extra=c(
    "Area",
    "Gap Fill Status",
    "Peak Rating"),
 mz_colname="m/z",
 rt_colname="RT [min]",
 type=c("pos",
    "neg"),
 se_list=NULL,
 ...)
{
   type <- match.arg(type);

   # extract annotation and assay data
   names(assay_names_extra) <- gsub("[ ]+",
      "_",
      tolower(assay_names_extra))
   non_assay_colnames <- jamba::unvigrep(
      paste0("^",
         jamba::cPaste(assay_names_extra, "|"),
         ":"),
      colnames(x))
   # head(x[, non_assay_colnames], 5)

   x_names <- paste0(
      seq_len(nrow(x)),
      "|",
      format(
         round(x[, mz_colname]*(10^nsmall_mz))/(10^nsmall_mz),
         nsmall=nsmall_mz,
         trim=TRUE),
      "|",
      format(
         round(x[, rt_colname]*(10^nsmall_rt))/(10^nsmall_rt),
         nsmall=nsmall_rt,
         trim=TRUE))
   # print(jamba::vigrep("^(27|12)[|]", x_names));
   rownames(x) <- x_names;

   # prepare rowData
   x_rowData <- x[, non_assay_colnames, drop=FALSE];

   # assays
   assays_list <- lapply(assay_names_extra, function(iassay_name){
      iassay_colnames <- jamba::vigrep(
         paste0("^", iassay_name, ":"), colnames(x))
      if (length(iassay_colnames) == 0) {
         return(NULL)
      }
      ix <- x[, iassay_colnames, drop=FALSE]
      # remove prefix
      colnames(ix) <- gsub(paste0("^", iassay_name, ":[ ]*"),
         "",
         colnames(ix))
      # remove suffix
      colnames(ix) <- gsub("[ ]*[(][A-Z0-9]+[)][ ]*$",
         "",
         colnames(ix))
      # remove suffix
      colnames(ix) <- gsub("[.](raw)$",
         "",
         colnames(ix))
      ix <- as.matrix(ix)
      if (max(ix, na.rm=TRUE) > 1000) {
         ix <- log2(1 + ix);
      }
      ix
   })
   # remove empty entries
   assays_list <- jamba::rmNULL(assays_list);
   if (length(assays_list) == 0) {
      jamba::printDebug("process_metab_compounds_file(): ",
         "No matching assay names were recognized: ",
         paste0("'", assay_names_extra, "'"),
         ". Returning NULL.")
      return(se_list)
   }
   se <- SummarizedExperiment::SummarizedExperiment(
      assays=assays_list,
      rowData=x_rowData)

   # now merge with incoming se
   # match rows to se_list
   if (length(se_list) > 0) {
      se_names <- jamba::vigrep(type, names(se_list));
      for (se_name in se_names) {
         # match by rownum only - no other matching is reliable
         match_pos_rows <- match(
            gsub("[|].+$", "",
               rownames(se_list[[se_name]])),
            gsub("[|].+$", "",
               rownames(se)))
         rownames(se)[match_pos_rows] <- rownames(se_list[[se_name]])
         match_pos_columns <- match(colnames(se_list[[se_name]]),
            colnames(se))
         for (iassay_name in SummarizedExperiment::assayNames(se)) {
            SummarizedExperiment::assays(
               se_list[[se_name]])[[iassay_name]] <- (
                  SummarizedExperiment::assays(
                     se[match_pos_rows, match_pos_columns])[[iassay_name]])
         }
      }
      return(se_list)
   }

   return(se)
}


#' Convert imputed metabolomics assay data to NA
#'
#' Convert imputed metabolomics assay data to NA
#'
#' This function takes `SummarizedExperiment` data as input, which
#' must contain at least two `SummarizedExperiment::assayNames()`
#' entries:
#' 1. `numeric` assay measurements
#' 2. `numeric` impute flag data
#'
#' The coordinates of cells in the impute flag data
#' which match `impute_flags` are used to convert data in assay measurements
#' to `NA`.
#'
#' @family jam utility functions
#'
#' @returns `SummarizedExperiment` or when input `se` is a `list`
#'    a `list` of `SummarizedExperiment` objects is returned.
#'
#' @param se `SummarizedExperiment` object, or `list` of
#'    `SummarizedExperiment` objects.
#'    When a `list` is provided, each entry is iterated and processed,
#'    and a `list` is returned.
#' @param measurement_assay_names `character` vector with one or more
#'    `SummarizedExperiment::assayNames(se)` for which the imputed
#'    data will be converted to `NA` values.
#' @param impute_assay_name `character` string with one
#'    `SummarizedExperiment::assayNames(se)` entry that contains the
#'    impute flags. Only the first matching assay name is used.
#' @param suffix `character` suffix appended to each value in
#'    `measurement_assay_names` to create a new assay name for the
#'    resulting data matrix.
#'
#'    When filtering by different impute criteria, it may be helpful
#'    to include the codes in the suffix, for example filtering by only
#'    `impute_flags=8` may suggest `suffix="_noimpute8"`.
#' @param impute_flags `numeric` vector of flags indicating values in
#'    `impute_assay_name` assay data which represent imputed data.
#'
#'    The defaults:
#'    * 8: filled by trace area (imputed)
#'    * 32: filled by spectrum noise (imputed)
#'
#'    Other values:
#'    * 0: original ion used
#'    * 64: filled by matching ion
#'    * 128: filled by re-detected peak
#' @param ... additional arguments are ignored.
#'
#' @export
convert_imputed_assays_to_na <- function
(se,
 measurement_assay_names="area",
 impute_assay_name="gap_fill_status",
 suffix="_noimpute",
 impute_flags=c(8, 32),
 ...)
{
   #
   if (is.list(se)) {
      for (i in seq_along(se)) {
         se[[i]] <- convert_imputed_assays_to_na(se[[i]],
            measurement_assay_names=measurement_assay_names,
            impute_assay_name=impute_assay_name,
            suffix=suffix,
            impute_flags=impute_flags,
            ...)
      }
      return(se)
   }

   if (length(impute_flags) == 0) {
      stop("No impute_flags were defined.")
   }

   use_impute_assay_name <- head(intersect(impute_assay_name,
      SummarizedExperiment::assayNames(se)), 1)
   if (length(use_impute_assay_name) == 0) {
      stop(paste0("impute_assay_name '",
         impute_assay_name,
         "' was not found in assayNames(se)."))
   }

   use_assay_names <- intersect(measurement_assay_names,
      SummarizedExperiment::assayNames(se))
   if (length(use_assay_names) == 0) {
      stop(paste0("No measurement_assay_names (",
         jamba::cPaste(paste0("'", measurement_assay_names, "'")),
         ") were present in ",
         "assayNames(se)."))
   }

   # define imputed data positions
   is_imputed <- SummarizedExperiment::assays(se)[[use_impute_assay_name]]
   is_imputed[] <- (SummarizedExperiment::assays(se)[[use_impute_assay_name]]
      %in% c(impute_flags))

   for (use_assay_name in use_assay_names) {
      new_assay_name <- paste0(use_assay_name, suffix)
      SummarizedExperiment::assays(se)[[new_assay_name]] <- (
         SummarizedExperiment::assays(se)[[use_assay_name]]);
      if (any(is_imputed %in% 1)) {
         SummarizedExperiment::assays(se)[[
            new_assay_name]][is_imputed %in% 1] <- NA;
      }
   }
   return(se);
}
