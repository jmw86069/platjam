
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
      metapaths <- jamba::rmNA(get_salmon_root(metafile));
      metapaths <- metapaths[!duplicated(metapaths)];
      metajsons <- lapply(seq_along(metapaths), function(i){
         idf <- get_salmon_meta(metapaths[i],
            exclude_hashes=exclude_hashes,
            ...);
         idf$temp_rowname <- rownames(idf);
         idf;
      });
      metajson <- jamba::mergeAllXY(metajsons);
      rownames(metajson) <- metajson$temp_rowname;
      metajson <- metajson[, setdiff(colnames(metajson),
         "temp_rowname"), drop=FALSE];
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
   flenfile <- file.path(metapath,
      "libParams/flenDist.txt");

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
   rownames(metajson) <- names(metapath);
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

   # fragment length file
   if (file.exists(flenfile)) {
      flen_values <- parse_salmon_flenfile(flenfile,
         do_plot=FALSE);
      names(flen_values) <- paste0("fraglen_", names(flen_values))
      metajson[,names(flen_values)] <- flen_values;
   }

   return(metajson);
}

#' Parse Salmon fragment length file
#'
#' Parse Salmon fragment length file
#'
#' @param x `character` path to file usually named `"flenDist.txt"`.
#' @param k `numeric` used with `caTools::runmean()` to define the
#'    width for smoothing a running mean, helpful to de-noise
#'    the fragment length profile.
#' @param do_plot `logical` indicating whether to plot the data and
#'    resulting values.
#' @param retcolors `character` vector with three colors used when
#'    `do_plot=TRUE`: mean, median, mode.
#' @param ... additional arguments are passed to `plot()` when `do_plot=TRUE`.
#'
#' @family jam nextgen sequence functions
#'
#' @return `numeric` vector with summary values:
#' * mean: the weighted mean across the fragment length distribution
#' * median: the position with 50% cumulative density at or below
#' * mode: the position with the highest density
#' * each value in `probs` representing percentiles
#'
#' value representing the weighted mean across the
#'    fragment lengths reported, typically values from 1 to 1000.
#'
#' @export
parse_salmon_flenfile <- function
(x,
 k=7,
 do_plot=FALSE,
 retcolors=c("red",
    "gold",
    "dodgerblue",
    "grey"),
 probs=c(0.25, 0.75),
 ...)
{
   if (is.numeric(x)) {
      xvals <- x;
      xlens <- seq_along(x);
   } else if (!file.exists(x)) {
      stop(paste0("Salmon fragment length file is not accessible: ", x));
   } else {
      xdf <- data.table::fread(x,
         data.table=FALSE,
         header=FALSE);
      xvals <- unname(unlist(xdf[1,]));
      xlens <- seq_along(xvals)
   }
   names(xvals) <- xlens;
   xmean <- weighted.mean(x=xlens,
      w=xvals,
      na.rm=TRUE);
   if (jamba::check_pkg_installed("caTools")) {
      xsmooth <- caTools::runmean(xvals,
         k=k,
         endrule="mean");
      xmode <- xlens[which.max(xsmooth)];
   } else {
      xmode <- xlens[which.max(xvals)];
      xsmooth <- NULL;
   }
   w <- xvals * (1 / sum(xvals));
   probs <- setdiff(probs, 0.5);
   xmedian <- xlens[match(TRUE, cumsum(w) >= 0.5)]
   xprobs <- sapply(probs, function(i){
      xlens[match(TRUE, cumsum(w) >= i)]
   })
   names(xprobs) <- paste0(round(100 * probs), "%");

   retvals <- c(mean=xmean,
      median=xmedian,
      mode=xmode,
      xprobs);

   # optional plot
   if (do_plot) {
      plot(x=xlens,
         y=xvals,
         type="l",
         col="blue",
         lwd=4,
         ...);
      if (length(xsmooth) > 0) {
         lines(x=xlens,
            y=xsmooth,
            col="red",
            lwd=2);
      }
      abline(v=retvals,
         lwd=1,
         col=rep(retcolors,
            c(1, 1, 1, length(xprobs))),
         lty="dashed")
      legend("topright",
         legend=names(retvals),
         col=rep(retcolors,
            c(1, 1, 1, length(xprobs))),
         lwd=2,
         lty=1)
   }
   retvals
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
#' @param ... additional arguments are ignored.
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
(file,
 ...)
{
   ## Try to use rprojroot to find the Salmon root directory
   if (!suppressPackageStartupMessages(require(rprojroot))) {
      stop("The rprojroot package is required for get_salmon_meta().");
   }
   if (length(file) > 1) {
      salmonroots <- sapply(file, get_salmon_root);
      return(salmonroots);
   }
   names_file <- names(file);
   salmonroot <- tryCatch({
      rprojroot::find_root(path=file,
         "cmd_info.json");
   }, error=function(e){
      NA;
   });
   names(salmonroot) <- names_file;
   return(salmonroot);
}

#' Save Salmon QC to Xlsx
#'
#' @family jam nextgen sequence functions
#'
#' @param metajsons `data.frame` output from `get_salmon_meta()`
#' @param salmon_qc_xlsx `character` path to file to be saved.
#' @param adj_target_reads `numeric` number indicating the target
#'    number of sequence reads per row, used to define color gradient
#'    range.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to underlying functions.
#'
#' @export
save_salmon_qc_xlsx <- function
(metajsons,
 salmon_qc_xlsx=NULL,
 adj_target_reads=75000000,
 verbose=FALSE,
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
   qc_cols <- jamba::provigrep(qc_cols,
      colnames(metajsons));
   meta_qc <- metajsons[,qc_cols,drop=FALSE];
   colorSub1 <- jamba::nameVector(
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
      jamba::writeOpenxlsx(file=salmon_qc_xlsx,
         x=meta_qc,
         sheetName="Raw_QC",
         colorSub=colorSub1,
         intColumns=unname(which(colMeans(meta_qc[,-1]) >= 10) + 1),
         numColumns=unname(which(colMeans(meta_qc[,-1]) < 10) + 1),
         numFormat="#,##0.00000",
         verbose=verbose,
         doConditional=FALSE);
      applyXlsxConditionalFormatByColumn(file=salmon_qc_xlsx,
         x=meta_qc,
         sheet="Raw_QC",
         dryrun=FALSE,
         verbose=verbose,
         numColumns=seq(from=2, to=ncol(meta_qc))
      )
      jamba::set_xlsx_colwidths(salmon_qc_xlsx,
         sheet="Raw_QC",
         widths=rep(c(50,12),c(1, ncol(meta_qc)-1)))
      jamba::set_xlsx_rowheights(salmon_qc_xlsx,
         sheet="Raw_QC",
         rows=seq_len(nrow(meta_qc)+1),
         heights=rep(c(17*5,17),c(1, nrow(meta_qc))))

      ## adjusted QC
      jamba::writeOpenxlsx(file=salmon_qc_xlsx,
         x=metam_adj_df,
         sheetName="Adjusted_QC",
         colorSub=colorSub1,
         intColumns=unname(which(colMeans(metam_adj_df[,-1]) >= 10) + 1),
         numColumns=unname(which(colMeans(metam_adj_df[,-1]) < 10) + 1),
         numFormat="#,##0.00",
         append=TRUE,
         verbose=verbose,
         doConditional=FALSE);
      applyXlsxConditionalFormatByColumn(file=salmon_qc_xlsx,
         x=metam_adj_df,
         sheet="Adjusted_QC",
         dryrun=FALSE,
         verbose=verbose,
         numColumns=seq(from=2, to=ncol(metam_adj_df))
      )
      jamba::set_xlsx_colwidths(salmon_qc_xlsx,
         sheet="Adjusted_QC",
         widths=rep(c(50,12),c(1, ncol(metam_adj_df)-1)))
      jamba::set_xlsx_rowheights(salmon_qc_xlsx,
         sheet="Adjusted_QC",
         rows=seq_len(nrow(metam_adj_df)+1),
         heights=rep(c(17*5,17),c(1, nrow(metam_adj_df))))

      ## Salmon orientation
      ori_means <- colMeans(meta_ori_pct_df[,-1,drop=FALSE]);
      ori_int <- (ori_means == 0 | ori_means >= 10);
      jamba::writeOpenxlsx(file=salmon_qc_xlsx,
         x=meta_ori_df,
         sheetName="Orientation",
         colorSub=colorSub1,
         intColumns=unname(which(ori_int) + 1),
         numColumns=unname(which(!ori_int) + 1),
         numFormat="#,##0.0",
         append=TRUE,
         verbose=verbose,
         doConditional=FALSE);
      applyXlsxConditionalFormatByColumn(file=salmon_qc_xlsx,
         x=meta_ori_df,
         sheet="Orientation",
         dryrun=FALSE,
         verbose=verbose,
         numColumns=seq(from=2, to=ncol(meta_ori_df))
      )
      jamba::set_xlsx_colwidths(salmon_qc_xlsx,
         sheet="Orientation",
         widths=rep(c(50,12),c(1, ncol(meta_ori_df)-1)))
      jamba::set_xlsx_rowheights(salmon_qc_xlsx,
         sheet="Orientation",
         rows=seq_len(nrow(meta_ori_df)+1),
         heights=rep(c(17*5,17),c(1, nrow(meta_ori_df))))

      ## Salmon orientation by percent
      ori_pct_means <- colMeans(meta_ori_pct_df[,-1]);
      ori_pct_int <- (ori_pct_means == 0 | ori_pct_means >= 10);
      jamba::writeOpenxlsx(file=salmon_qc_xlsx,
         x=meta_ori_pct_df,
         sheetName="Orientation_Percent",
         colorSub=colorSub1,
         intColumns=unname(which(ori_pct_int) + 1),
         numColumns=unname(which(!ori_pct_int) + 1),
         numFormat="#,##0.000",
         append=TRUE,
         verbose=verbose,
         doConditional=FALSE);
      applyXlsxConditionalFormatByColumn(file=salmon_qc_xlsx,
         x=meta_ori_pct_df,
         sheet="Orientation_Percent",
         dryrun=FALSE,
         verbose=verbose,
         numColumns=seq(from=2, to=ncol(meta_ori_pct_df))
      )
      jamba::set_xlsx_colwidths(salmon_qc_xlsx,
         sheet="Orientation_Percent",
         widths=rep(c(50,12),c(1, ncol(meta_ori_pct_df)-1)))
      jamba::set_xlsx_rowheights(salmon_qc_xlsx,
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
         jamba::printDebug("applyXlsxConditionalFormatByColumn(): ",
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
         jamba::applyXlsxConditionalFormat(xlsxFile=file,
            sheet=sheet,
            numStyle=numStyle,
            numRule=x_rule,
            numColumns=i,
            ...);
      }
   }
}

