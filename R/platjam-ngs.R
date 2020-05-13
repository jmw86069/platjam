
#' Get Salmon metadata and aux info into a data.frame
#'
#' Get Salmon metadata and aux info into a data.frame
#'
#' This function takes a file path to one or more Salmon
#' output files, uses that path to locate the full set of
#' available files, loads data from each of the discovered
#' files, and returns the results in a `data.frame` format.
#'
#' This function uses `rprojroot::find_root()` to find the
#' root directory, defined as the directory that contains
#' the file `"cmd_info.json"`. The path to `"meta_info.json"`
#' is constructed relative to that location.
#'
#' Recognized files:
#'
#' * `meta_info.json` - typically in a subdirectory `aux_info/meta_info.json`
#' * `cmd_info.json` - typically in the same directory as the `aux_info`
#' directory.
#'
#' If a relative path to `"cmd_info.json"` cannot be determined, this
#' function returns `NULL`.
#'
#' When the input `metafile` includes multiple files, only
#' the unique Salmon root directories are returned.
#'
#' This function uses `jsonlite` to read each JSON file, which
#' is converted to a `data.frame`. Any JSON fields that contain
#' multiple values are comma-delimited using `jamba::cPaste()`
#' in order to fit on one row in the `data.frame`.
#'
#' @family jam nextgen sequence functions
#'
#' @return `data.frame` whose number of rows is equal to the number
#'    of unique Salmon root directories in the input `metafile`.
#'    For any input `metafile` not found, the output is `NULL`.
#'
#' @param metafile character vector of one or more files, usually the
#'    full file path to the `meta_info.json` file after running Salmon
#'    quant. The path `metafile` should be the path to any output file
#'    from one Salmon quant analysis.
#' @param exclude_hashes logical indicating whether to drop columns that
#'    contain file hashes.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' cmdinfopath <- system.file("data", "salmonOut", "cmd_info.json", package="platjam");
#' if (nchar(cmdinfopath) > 0) {
#'    get_salmon_meta(cmdinfopath);
#' }
#'
#' @export
get_salmon_meta <- function
(metafile,
 exclude_hashes=TRUE,
 ...)
{
   ##
   if (!suppressPackageStartupMessages(require(rprojroot))) {
      stop("The rprojroot package is required for get_salmon_meta().");
   }
   if (!suppressPackageStartupMessages(require(jsonlite))) {
      stop("The jsonlite package is required for get_salmon_meta().");
   }
   if (length(metafile) > 1) {
      metapaths <- unique(jamba::rmNA(get_salmon_root(metafile)));
      metajsons <- lapply(seq_along(metapaths), function(i){
         get_salmon_meta(metapaths[i],
            exclude_hashes=exclude_hashes,
            ...);
      });
      metajson <- jamba::mergeAllXY(metajsons);
      return(metajson);
   }

   ## Try to use rprojroot to find the Salmon root directory
   metapath <- get_salmon_root(metafile);
   if (is.na(metapath)) {
      return(NULL);
   }
   metafile <- file.path(metapath, "aux_info",
      "meta_info.json");
   cmdinfofile <- file.path(metapath,
      "cmd_info.json");
   libfile <- file.path(metapath,
      "lib_format_counts.json");

   json1 <- jsonlite::read_json(metafile);
   jsoninfo <- jamba::sdim(json1);
   if (any(jsoninfo$class %in% "list" & jsoninfo$rows > 0)) {
      jsonlist <- which(jsoninfo$class %in% "list" & jsoninfo$rows > 0);
      json1[jsonlist] <- jamba::cPaste(json1[jsonlist]);
   }
   if (length(names(metafile)) > 0) {
      metajson <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         name=names(metafile),
         as.data.frame(json1[jsoninfo$rows > 0]));
   } else {
      metajson <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         as.data.frame(json1[jsoninfo$rows > 0]));
   }
   if (exclude_hashes) {
      metajson <- metajson[,jamba::unvigrep("hash", colnames(metajson)),drop=FALSE];
   }

   ## cmd_info
   if (file.exists(cmdinfofile)) {
      json1 <- jsonlite::read_json(cmdinfofile);
      jsoninfo <- jamba::sdim(json1);
      if (any(jsoninfo$class %in% "list" & jsoninfo$rows > 0)) {
         jsonlist <- which(jsoninfo$class %in% "list" & jsoninfo$rows > 0);
         json1[jsonlist] <- jamba::cPaste(json1[jsonlist]);
      }
      if (any(jsoninfo$class %in% "list" & jsoninfo$rows == 0)) {
         jsonlist0 <- which(jsoninfo$class %in% "list" & jsoninfo$rows == 0);
         json1[jsonlist0] <- sapply(names(json1[jsonlist0]), function(i){i});
      }
      cmdjson <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         as.data.frame(json1));
      metajson[,colnames(cmdjson)] <- cmdjson;
   }

   ## lib_format_counts
   if (file.exists(libfile)) {
      json1 <- jsonlite::read_json(libfile);
      jsoninfo <- jamba::sdim(json1);
      if (any(jsoninfo$class %in% "list" & jsoninfo$rows > 0)) {
         jsonlist <- which(jsoninfo$class %in% "list" & jsoninfo$rows > 0);
         json1[jsonlist] <- jamba::cPaste(json1[jsonlist]);
      }
      if (any(jsoninfo$class %in% "list" & jsoninfo$rows == 0)) {
         jsonlist0 <- which(jsoninfo$class %in% "list" & jsoninfo$rows == 0);
         json1[jsonlist0] <- sapply(names(json1[jsonlist0]), function(i){i});
      }
      libjson <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         as.data.frame(json1[jsoninfo$rows > 0]));
      metajson[,colnames(libjson)] <- libjson;
   }

   return(metajson);
}

#' Get Salmon root directory
#'
#' Get Salmon root directory
#'
#' This function uses the `rprojroot` package to find the
#' root output directory of a Salmon quant analysis, looking
#' for the file `"cmd_info.json"`. If this file is not found,
#' this function returns `NULL`.
#'
#' @return file path to the Salmon root directory, or `NA` if
#'    the Salmon root directory could not be found.
#'
#' @param file character path to a Salmon output file, or a vector
#'    of files.
#'
#' @family jam nextgen sequence functions
#'
#' @examples
#' cmdinfopath <- system.file("data", "salmonOut", "cmd_info.json", package="platjam");
#' if (nchar(cmdinfopath) > 0) {
#'    get_salmon_root(cmdinfopath);
#' }
#'
#' @export
get_salmon_root <- function
(file)
{
   ## Try to use rprojroot to find the Salmon root directory
   if (!suppressPackageStartupMessages(require(rprojroot))) {
      stop("The rprojroot package is required for get_salmon_meta().");
   }
   if (length(file) > 1) {
      salmonroots <- sapply(file, get_salmon_root);
      return(salmonroots);
   }
   salmonroot <- tryCatch({
      rprojroot::find_root(path=file,
         "cmd_info.json");
   }, error=function(e){
      NA;
   });
   return(salmonroot);
}

#' Save Salmon QC to Xlsx
#'
#' @family jam nextgen sequence functions
#'
#' @export
save_salmon_qc_xlsx <- function
(metajsons,
 salmon_qc_xlsx=NULL,
 adj_target_reads=75000000,
 ...)
{
   rownames(metajsons) <- metajsons$output;

   qc_cols <- c("output",
      "num_eq_classes",
      "num_processed",
      "num_mapped",
      "num_decoy_fragments",
      "num_dovetail_fragments",
      "num_fragments_filtered_vm",
      "num_alignments_below_threshold_for_mapped_fragments_vm",
      "percent_mapped",
      "num_compatible_fragments","num_assigned_fragments",
      "num_frags_with_concordant_consistent_mappings",
      "num_frags_with_inconsistent_or_orphan_mappings",
      "strand_mapping_bias");
   qc_cols <- provigrep(qc_cols,
      colnames(metajsons));
   meta_qc <- metajsons[,qc_cols,drop=FALSE];
   colorSub1 <- nameVector(
      colorjam::group2colors(gsub("[.].+", "", meta_qc[,1])),
      meta_qc[,1]);

   ## adjust to 75 million input reads
   qc_cols_adj <- unvigrep("output|bias", qc_cols);
   metam <- as.matrix(metajsons[,qc_cols_adj,drop=FALSE]);
   metam_adj <- metam;
   metam_adjustment <- adj_target_reads / metam[,"num_processed"];
   adj_cols <- jamba::unvigrep("percent", colnames(metam));
   metam_adj[,adj_cols] <- metam[,adj_cols,drop=FALSE] * metam_adjustment;
   metam_adj <- cbind(metam_adj, adjustment=metam_adjustment);
   metam_adj_df <- data.frame(check.names=FALSE,
      output=rownames(metam_adj),
      metam_adj);

   ## orientation data
   ori_cols <- jamba::vgrep("^[A-Z]+$", colnames(metajsons));
   meta_ori <- as.matrix(metajsons[,ori_cols,drop=FALSE]);
   meta_ori_df <- data.frame(check.names=FALSE,
      output=rownames(meta_ori),
      meta_ori);

   ## orientation data by percent for each sample
   meta_ori_pct <- meta_ori / rowSums(meta_ori) * 100;
   meta_ori_pct_df <- data.frame(check.names=FALSE,
      output=rownames(meta_ori_pct),
      meta_ori_pct);

   if (length(salmon_qc_xlsx) > 0) {
      writeOpenxlsx(file=salmon_qc_xlsx,
         x=meta_qc,
         sheetName="Raw_QC",
         colorSub=colorSub1,
         intColumns=unname(which(colMeans(meta_qc[,-1]) >= 10) + 1),
         numColumns=unname(which(colMeans(meta_qc[,-1]) < 10) + 1),
         numFormat="#,##0.00000",
         doConditional=FALSE);
      applyXlsxConditionalFormatByColumn(file=salmon_qc_xlsx,
         x=meta_qc,
         sheet="Raw_QC",
         dryrun=FALSE,
         numColumns=seq(from=2, to=ncol(meta_qc))
      )
      set_xlsx_colwidths(salmon_qc_xlsx,
         sheet="Raw_QC",
         widths=rep(c(50,12),c(1, ncol(meta_qc)-1)))
      set_xlsx_rowheights(salmon_qc_xlsx,
         sheet="Raw_QC",
         rows=seq_len(nrow(meta_qc)+1),
         heights=rep(c(17*5,17),c(1, nrow(meta_qc))))

      ## adjusted QC
      writeOpenxlsx(file=salmon_qc_xlsx,
         x=metam_adj_df,
         sheetName="Adjusted_QC",
         colorSub=colorSub1,
         intColumns=unname(which(colMeans(metam_adj_df[,-1]) >= 10) + 1),
         numColumns=unname(which(colMeans(metam_adj_df[,-1]) < 10) + 1),
         numFormat="#,##0.00",
         append=TRUE,
         doConditional=FALSE);
      applyXlsxConditionalFormatByColumn(file=salmon_qc_xlsx,
         x=metam_adj_df,
         sheet="Adjusted_QC",
         dryrun=FALSE,
         verbose=TRUE,
         numColumns=seq(from=2, to=ncol(metam_adj_df))
      )
      set_xlsx_colwidths(salmon_qc_xlsx,
         sheet="Adjusted_QC",
         widths=rep(c(50,12),c(1, ncol(metam_adj_df)-1)))
      set_xlsx_rowheights(salmon_qc_xlsx,
         sheet="Adjusted_QC",
         rows=seq_len(nrow(metam_adj_df)+1),
         heights=rep(c(17*5,17),c(1, nrow(metam_adj_df))))

      ## Salmon orientation
      ori_means <- colMeans(meta_ori_pct_df[,-1,drop=FALSE]);
      ori_int <- (ori_means == 0 | ori_means >= 10);
      writeOpenxlsx(file=salmon_qc_xlsx,
         x=meta_ori_df,
         sheetName="Orientation",
         colorSub=colorSub1,
         intColumns=unname(which(ori_int) + 1),
         numColumns=unname(which(!ori_int) + 1),
         numFormat="#,##0.0",
         append=TRUE,
         doConditional=FALSE);
      applyXlsxConditionalFormatByColumn(file=salmon_qc_xlsx,
         x=meta_ori_df,
         sheet="Orientation",
         dryrun=FALSE,
         verbose=TRUE,
         numColumns=seq(from=2, to=ncol(meta_ori_df))
      )
      set_xlsx_colwidths(salmon_qc_xlsx,
         sheet="Orientation",
         widths=rep(c(50,12),c(1, ncol(meta_ori_df)-1)))
      set_xlsx_rowheights(salmon_qc_xlsx,
         sheet="Orientation",
         rows=seq_len(nrow(meta_ori_df)+1),
         heights=rep(c(17*5,17),c(1, nrow(meta_ori_df))))

      ## Salmon orientation by percent
      ori_pct_means <- colMeans(meta_ori_pct_df[,-1]);
      ori_pct_int <- (ori_pct_means == 0 | ori_pct_means >= 10);
      writeOpenxlsx(file=salmon_qc_xlsx,
         x=meta_ori_pct_df,
         sheetName="Orientation_Percent",
         colorSub=colorSub1,
         intColumns=unname(which(ori_pct_int) + 1),
         numColumns=unname(which(!ori_pct_int) + 1),
         numFormat="#,##0.000",
         append=TRUE,
         doConditional=FALSE);
      applyXlsxConditionalFormatByColumn(file=salmon_qc_xlsx,
         x=meta_ori_pct_df,
         sheet="Orientation_Percent",
         dryrun=FALSE,
         verbose=TRUE,
         numColumns=seq(from=2, to=ncol(meta_ori_pct_df))
      )
      set_xlsx_colwidths(salmon_qc_xlsx,
         sheet="Orientation_Percent",
         widths=rep(c(50,12),c(1, ncol(meta_ori_pct_df)-1)))
      set_xlsx_rowheights(salmon_qc_xlsx,
         sheet="Orientation_Percent",
         rows=seq_len(nrow(meta_ori_pct_df)+1),
         heights=rep(c(17*5,17),c(1, nrow(meta_ori_pct_df))))
      dfs <- list(meta_qc=meta_qc,
         meta_qc_adj=metam_adj_df,
         meta_ori=meta_ori_df,
         meta_ori_pct=meta_ori_pct_df);
   }

   invisible(dfs);
}

#' Apply Xlsx Conditional Formatting by Column
#'
#' Apply Xlsx Conditional Formatting by Column
#'
#' This function applies conditional formatting to `.xlsx`
#' files, using the numeric range of each column to define
#' the color rules. It is an extension to
#' `jamba::applyXlsxConditionalFormat()` that differs only
#' in applying colors to the data range, rather than using
#' a pre-defined fixed range across all columns.
#'
#' @family jam export functions
#'
#' @export
applyXlsxConditionalFormatByColumn <- function
(file,
 x,
 sheet=1,
 numColumns=NULL,
 numStyle=c("#F2F2F2","#6464FC","#2929D4"),
 mid_rule=c("mid","mean","median"),
 dryrun=FALSE,
 verbose=FALSE,
 ...)
{
   mid_rule <- match.arg(mid_rule);
   for (i in numColumns) {
      y <- x[,i];
      x_range <- range(y, na.rm=TRUE);
      if (diff(x_range) == 0) {
         x_range <- range(c(0, 1, x_range));
      }
      if ("mid" %in% mid_rule) {
         x_mid <- mean(x_range);
      } else if ("mean" %in% mid_rule) {
         x_mid <- mean(y, na.rm=TRUE);
      } else if ("median" %in% mid_rule) {
         x_mid <- median(y, na.rm=TRUE);
      }
      x_rule <- c(x_range, x_mid)[c(1,3,2)];
      if (verbose || dryrun) {
         printDebug("applyXlsxConditionalFormatByColumn(): ",
            paste0(colnames(x)[i],": "),
            format(x_rule,
               digits=1,
               big.mark=",",
               scientific=FALSE,
               trim=TRUE),
            fgText=list("darkorange2",
               "dodgerblue",
               setTextContrastColor(numStyle)),
            bgText=list(NA, NA, numStyle));
      }
      if (!dryrun) {
         applyXlsxConditionalFormat(xlsxFile=file,
            sheet=sheet,
            numStyle=numStyle,
            numRule=x_rule,
            numColumns=i,
            ...);
      }
   }
}

#' Set column widths in Xlsx files
#'
#' Set column widths in Xlsx files
#'
#' This function is a light wrapper to perform these steps
#' from the very useful `openxlsx` R package:
#'
#' * `openxlsx::loadWorkbook()`
#' * `openxlsx::setColWidths()`
#' * `openxlsx::saveWorkbook()`
#'
#' @family jam export functions
#'
#' @export
set_xlsx_colwidths <- function
(xlsxFile,
   sheet=1,
   cols=seq_along(widths),
   widths=11,
   ...)
{
   ## Load the requested file as a workbook
   wb <- openxlsx::loadWorkbook(xlsxFile);

   openxlsx::setColWidths(wb,
      sheet=sheet,
      cols=cols,
      widths=widths,
      ...);

   ## Save workbook
   openxlsx::saveWorkbook(wb,
      xlsxFile,
      overwrite=TRUE);
}

#' Set row heights in Xlsx files
#'
#' This function is a light wrapper to perform these steps
#' from the very useful `openxlsx` R package:
#'
#' * `openxlsx::loadWorkbook()`
#' * `openxlsx::setRowHeights()`
#' * `openxlsx::saveWorkbook()`
#'
#' @family jam export functions
#'
#' @export
set_xlsx_rowheights <- function
(xlsxFile,
   sheet=1,
   rows=seq_along(heights)+1,
   heights=17,
   ...)
{
   ## Load the requested file as a workbook
   wb <- openxlsx::loadWorkbook(xlsxFile);

   openxlsx::setRowHeights(wb,
      sheet=sheet,
      rows=rows,
      heights=heights);

   ## Save workbook
   openxlsx::saveWorkbook(wb,
      xlsxFile,
      overwrite=TRUE);
}
