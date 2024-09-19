
#' Re-center coverage matrix data
#'
#' Re-center coverage matrix data
#'
#' Coverage matrix data is provided as `nmatlist` which is a
#' `list` of `normalizedMatrix` objects (see `EnrichmentHeatmap`).
#' One or more `recenter_heatmap` are defined, and the summit
#' is calculated for each row using `smooth.spline()`.
#'
#' For each row, the summit position is therefore interpolated.
#'
#' Coverage matrix data `nmatlist` is then shifted to recenter peaks
#' across all coverage files, and the summit offset is stored as
#' an attribute `attr(nmatlist, "summit_offset")`.
#'
#' ## Use Case
#'
#' The input represents sequence coverage data for a set of
#' genome regions of interest, for example ChIP-seq peaks, ATAC-seq peaks,
#' where it may be useful for the center of enrichment to
#' represent the position with highest coverage.
#' Many peak/regional enrichment calling tools may not provide
#' a summit position, the summit position may not be accurate,
#' and/or the summit across multiple sets of merged peaks may
#' not be available.
#'
#' This method can use coverage data across multiple `nmatlist`
#' matrix data to calculate a collective summit position.
#'
#' ## Recommendation
#'
#' The recommended workflow is to create coverage matrix data
#' for a region wider than used for the final figure, so that
#' the re-centering can be performed while maintaining coverage
#' throughout the desired range.
#'
#' @family jam coverage heatmap functions
#'
#' @returns object in format `nmatlist`, a `list` of `normalizedMatrix`
#'    objects.
#'
#' @param nmatlist `list` of `normalizedMatrix` objects
#' @param recenter_heatmap `numeric` (default 1) index with one or
#'    more entries in `nmatlist` to use for re-centering.
#' @param recenter_range `numeric` (default NULL) with optional
#'    maximum distance from the target (center) of coverage
#'    in `nmatlist`. For example, if `nmatlist` data spans -50kb to +50kb,
#'    but peaks are no wider than 1kb,  consider using
#'    `recenter_range=1000` so that the recentering will only use
#'    coverage data -1000bp to +1000bp at most.
#' @param recenter_invert `logical` indicating whether to invert the
#'    coverage, therefore effectively taking the minimum signal.
#'    This value is recycled to `length(recenter_heatmap)` such that
#'    each heatmap can individually be inverted as relevant.
#' @param spar.edge.buffer `numeric` values passed to `summit_from_vector()`
#' @param empty_value `numeric` value used for empty values created
#'    by the "edges" of recentered matrix data. Default is 0, other
#'    values may not be well-supported.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' ## There is a small example file to use for testing
#' # library(jamba)
#' cov_file1 <- system.file("data", "tss_coverage.matrix", package="platjam");
#' cov_file2 <- system.file("data", "h3k4me1_coverage.matrix", package="platjam");
#' cov_files <- c(cov_file1, cov_file2);
#' names(cov_files) <- gsub("[.]matrix",
#'    "",
#'    basename(cov_files));
#' nmatlist <- coverage_matrix2nmat(cov_files, verbose=FALSE);
#'
#' nmatlist2heatmaps(nmatlist,
#'    title="Input data",
#'    transform=c("log2signed", "sqrt"));
#'
#' nmatlist1 <- recenter_nmatlist(nmatlist)
#' nmatlist2heatmaps(nmatlist1,
#'    title="Input data, recentered by tss signal",
#'    transform=c("log2signed", "sqrt"));
#'
#' nmatlist2i <- recenter_nmatlist(nmatlist, recenter_heatmap=2, recenter_invert=TRUE)
#' nmatlist2heatmaps(nmatlist2i,
#'    title="Input data, recentered by inverted h3k4me1 signal",
#'    transform=c("log2signed", "sqrt"));
#' head(data.frame(summit_name=attr(nmatlist2i[[1]], "summit_name")))
#'
#' nmatlist2is <- restrand_nmatlist(nmatlist2i, restrand_heatmap=2, recenter_invert=FALSE)
#' nmatlist2heatmaps(nmatlist2is,
#'    title="Input data, recentered by inverted h3k4me1 signal,\nrestranded by tss",
#'    transform=c("log2signed", "sqrt"));
#'
#' # summarize recenter and restrand output
#' head(data.frame(
#'    row=attr(nmatlist2is[[1]], "dimnames")[[1]],
#'    summit_name=attr(nmatlist2is[[1]], "summit_name"),
#'    restrand=attr(nmatlist2is[[1]], "restrand")))
#'
#' @export
recenter_nmatlist <- function
(nmatlist,
 recenter_heatmap=1,
 recenter_range=NULL,
 recenter_invert=FALSE,
 spar=0.5,
 edge_buffer=0,
 empty_value=0,
 verbose=FALSE,
 ...)
{
   ## Basic steps:
   # - iterate each row in `nmatlist`
   # - define aggregate numeric vector across all `recenter_heatmap`
   # - summit_from_vector() and store new summit position
   # - re-center nmatlist using the new summit position

   ## Prep work
   #
   # Confirm all bin_width are identical
   # recenter_binwidths <- lapply(nmatlist[recenter_heatmap], function(nmat){
   #    if (!"extend" %in% names(attributes(nmat)) ||
   #          !"upstream_index" %in% names(attributes(nmat))) {
   #       # Note: if binwidth does not exist, it could be because this is
   #       # gene-body binned data
   #    }
   #    bin_width <- round(attr(nmat, "extend")[1] /  length(attr(nmat, "upstream_index")))
   # })

   ## Define the nmatlist with appropriate range to use for recentering
   if (length(recenter_invert) == 0) {
      recenter_invert <- FALSE;
   }
   recenter_invert <- rep(recenter_invert,
      length.out=length(recenter_heatmap));
   recenter_nmatlist <- lapply(nmatlist[recenter_heatmap], function(nmat){
      if (length(recenter_range) == 1 && recenter_range > 0) {
         nmat <- zoom_nmat(nmat,
            upstream_length=recenter_range,
            downstream_length=recenter_range,
            ...)
      }
      nmat;
   })
   nmatlist_nrow <- sapply(nmatlist, nrow);
   if (length(unique(nmatlist_nrow)) > 1) {
      stop("nrow must be identical across all nmatlist entries.")
   }

   ## Define which colnames to use
   recenter_ncols <- lapply(recenter_nmatlist, ncol)
   recenter_colnames_list <- lapply(recenter_nmatlist, colnames)
   # use the superset of colnames available
   recenter_colnames <- jamba::mixedSort(
      Reduce("union", recenter_colnames_list),
      keepNegative=TRUE);
   # expand colnames as needed
   if (length(unique(recenter_ncols)) > 1) {
      jamba::printDebug("recenter_nmatlist(): ",
         "Note that nmatlist contains inconsitent columns.",
         " Adjusting accordingly.");
   }
   # Actually, convert all to numeric matrix for convenience
   recenter_nmatlist <- lapply(seq_along(recenter_nmatlist), function(inum){
      nmat <- recenter_nmatlist[[inum]];
      match_col <- match(recenter_colnames, colnames(nmat));
      nmat_rownames <- rownames(nmat);
      nmat <- jamba::rmNA(naValue=0,
         nmat[, match_col, drop=FALSE]);
      colnames(nmat) <- recenter_colnames;
      rownames(nmat) <- nmat_rownames;
      if (TRUE %in% recenter_invert[[inum]]) {
         nmat <- max(nmat, na.rm=TRUE) - nmat;
      }
      nmat;
   });
   recenter_mat <- Reduce("+", recenter_nmatlist);

   ## Phase One: Determine summit position per row
   use_rows <- rownames(recenter_mat);
   new_summits <- data.frame(check.names=FALSE, jamba::rbindList(
      lapply(jamba::nameVector(use_rows), function(irow){
         x <- recenter_mat[irow, ]
         x_summit <- summit_from_vector(x,
            spar=spar,
            edge_buffer=edge_buffer,
            ...);
         # summit <- x_summit["summit"];
         # summit_height <- x_summit[["summit_height"]];
         x_summit;
      })))
   new_summits$summit_name <- recenter_colnames[new_summits$summit];
   ret_vals <- list();
   ret_vals$summits <- new_summits;


   ## Phase Two: Adjust each matrix
   new_nmatlist <- lapply(nmatlist, function(nmat){
      # even_ncol <- (ncol(nmat) %% 2) == 0;
      len1 <- ceiling(ncol(nmat) / 2);
      len2 <- floor(ncol(nmat) / 2);
      new_nmat <- jamba::rbindList(
         lapply(seq_along(new_summits$summit_name), function(inum){
            icol <- new_summits$summit_name[inum];
            match_icol <- match(icol, colnames(nmat));
            keep1 <- jamba::noiseFloor(
               seq(to=match_icol, by=1, length.out=len1),
               minimum=1,
               newValue=Inf)
            keep2 <- seq(from=match_icol + 1, by=1, length.out=len2);
            keep12 <- match(
               jamba::rmNA(
                  naValue=0,
                  colnames(nmat)[c(keep1, keep2)]),
               colnames(nmat))
            nmatline <- jamba::rmNA(
               naValue=empty_value,
               nmat[inum, keep12, drop=FALSE])
            colnames(nmatline) <- colnames(nmat)
            nmatline
         }))
      attributes(new_nmat) <- attributes(nmat)
      attr(new_nmat, "summit_name") <- new_summits$summit_name;
      new_nmat;
   })
   new_nmatlist

   # ret_vals
}



#' Determine summit from numeric vector
#'
#' Determine summit from numeric vector
#'
#' This function takes a numeric vector, intended to be data that
#' represents some signal across a range where that signal is above
#' noise; it calls `smooth.spline()` to generate a smooth curve across
#' the region, then returns the x position with the max smoothed
#' spline signal.
#'
#' The original intent is to take genome sequence coverage across
#' an enriched region (a "peak") and determine the peak summit.
#' It should work well for each row of a coverage matrix, provided the
#' coverage matrix is wide enough that the highest signal is located
#' inside the range analyzed.
#'
#' The other alternative is to import bigWig coverage data for a set
#' of regions of interest defined by a `GRanges` object.
#' A useful function is `splicejam::getGRcoverageFromBw()` which
#' can load coverage from one or multiple bigWig files, returning
#' a `GRanges` object with one column per bigWig file loaded.
#' Then iterate each coverage vector to determine the summit.
#'
#' @family jam utility functions
#'
#' @param x `numeric` vector from which a summit will be determined.
#' @param spar `numeric` or `NULL` passed to `smooth.spline()` to
#'    adjust the smoothing parameter. The default `spar=0.5` appears
#'    to provide smoothing at a reasonable and consistent level for
#'    genome coverage data, which tends to have long stretches of
#'    horizontal coverage that tend to be overfitted when
#'    `spar=NULL`.
#' @param edge_buffer `integer` number of values at the leading
#'    and trailing edge of `x` to be ignored when determining the
#'    summit. This argument is experimental, and is intended to
#'    prevent the very beginning or end of a region from being
#'    the "summit" when there may be an internal peak that is
#'    preferred. Note that when `(edge_buffer*2) > length(x)`
#'    the entire region is ignored, in which case the middle
#'    position is returned.
#' @param ... additional arguments are passed to `smooth.spline()`.
#'
#' @returns `integer` vector with two values:
#'    * `"summit"` with the index position of the highest point
#'    on the smoothed spline curve.
#'    If `x` has one uniform numeric value across the entire range,
#'    it returns the midpoint defined by `round(length(x)/2)`.
#'    If are two maximum values, the first position is returned.
#'    * `"summit_height"` `numeric` value with the spline height
#'    at the summit position.
#'
#' @export
summit_from_vector <- function
(x,
 spar=0.5,
 edge_buffer=0,
 return_height=TRUE,
 ...)
{
   # optional edge_buffer
   if (edge_buffer > 0) {
      minx <- min(x, na.rm=TRUE);
      if (length(x) > (edge_buffer * 2)) {
         xseq <- seq_along(x);
         x[head(xseq, edge_buffer)] <- minx;
         x[tail(xseq, edge_buffer)] <- minx;
      } else {
         x <- minx;
      }
   }
   # for x with no change, return the middle position
   if (length(unique(x)) == 1) {
      summitx <- round(length(x) / 2);
      summith <- head(x, 1);
   } else {
      # call smooth.spline
      smth <- smooth.spline(x=seq_along(x),
         y=x,
         spar=spar,
         ...);
      summitx <- which.max(smth$y);
      summith <- smth$y[summitx];
   }
   if (TRUE %in% return_height) {
      return(c(summit=summitx,
         summit_height=summith));
   } else {
      return(c(summit=summitx))
   }
}


#' Re-strand coverage matrix data by inferred strandedness
#'
#' Re-strand coverage matrix data by inferred strandedness
#'
#' Coverage matrix data is provided as `nmatlist` which is a
#' `list` of `normalizedMatrix` objects (see `EnrichmentHeatmap`).
#' One or more `restrand_heatmap` are defined, and the side with
#' higher coverage "wins" such that the downstream signal is
#' always higher within the specified range.
#'
#' @param nmatlist `list` of `normalizedMatrix` objects
#' @param restrand_heatmap `numeric` (default 1) index with one or
#'    more entries in `nmatlist` to use for re-stranding.
#' @param restrand_range `numeric` (default NULL) with optional
#'    maximum distance from the target (center) of coverage
#'    in `nmatlist`. For example, if `nmatlist` data spans -50kb to +50kb,
#'    but peaks are no wider than 1kb,  consider using
#'    `recenter_range=1000` so that the recentering will only use
#'    coverage data -1000bp to +1000bp at most.
#' @param restrand_buffer `numeric` not yet implemented, intended to
#'    enforce a "buffer" distance away from the center before
#'    calculating coverage.
#' @param restrand_invert `logical` indicating whether to invert the
#'    coverage, therefore effectively taking the minimum signal.
#'    This value is recycled to `length(restrand_heatmap)` such that
#'    each heatmap can individually be inverted as relevant.
#' @param empty_value `numeric` value used for empty values created
#'    by the "edges" of recentered matrix data. Default is 0, other
#'    values may not be well-supported.
#' @param transform `character` or `function` applied to `restrand_heatmap`
#'    to transform the data prior to calculating the coverage.
#'    * When using `transform` with `nmatlist2heatmaps()` it is recommended
#'    to use it here also.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @family jam coverage heatmap functions
#'
#' @returns object in format `nmatlist`, a `list` of `normalizedMatrix`
#'    objects.
#'
#' @examples
#' cov_file1 <- system.file("data", "tss_coverage.matrix", package="platjam");
#' cov_file2 <- system.file("data", "h3k4me1_coverage.matrix", package="platjam");
#' cov_files <- c(cov_file1, cov_file2);
#' names(cov_files) <- gsub("[.]matrix",
#'    "",
#'    basename(cov_files));
#' nmatlist <- coverage_matrix2nmat(cov_files, verbose=FALSE);
#'
#' # scramble the strandedness
#' new_k <- sample(c(TRUE, FALSE), size=nrow(nmatlist[[1]]), replace=TRUE)
#' nmatlist1 <- lapply(nmatlist, function(nmat){
#'    nmat[new_k, ] <- nmat[new_k, rev(colnames(nmat)), drop=FALSE]
#'    nmat;
#' })
#'
#' nmatlist2heatmaps(nmatlist1,
#'    title="Test data\n(mixed strandedness)",
#'    transform=c("log2signed", "sqrt"));
#'
#' nmatlist2heatmaps(nmatlist1,
#'    title="Test data,\nrestranded by TSS",
#'    restrand_heatmap=1,
#'    transform=c("log2signed", "sqrt"));
#'
#' nmatlist2heatmaps(nmatlist1,
#'    title="Test data,\nrecentered by H3K4me1",
#'    recenter_heatmap=2, recenter_invert=TRUE,
#'    transform=c("log2signed", "sqrt"));
#'
#' nmatlist2heatmaps(nmatlist1,
#'    title="Test data, recentered by H3K4me1,\nrestranded by TSS",
#'    recenter_heatmap=2, recenter_invert=TRUE,
#'    restrand_heatmap=1,
#'    transform=c("log2signed", "sqrt"));
#'
#' nmatlist_c <- restrand_nmatlist(nmatlist, restrand_invert=FALSE,
#'    transform=c("log2signed"), restrand_heatmap=1,
#'    restrand_range=600, restrand_buffer=100)
#' nmhm <- nmatlist2heatmaps(nmatlist_c,
#'    title="Re-stranded using tss signal",
#'    transform=c("log2signed", "sqrt"));
#'    k_clusters=8, k_method="euclidean", min_rows_per_k=10, k_heatmap=1,
#'
#' nmatlist2i <- recenter_nmatlist(nmatlist, recenter_heatmap=1, recenter_range=1000, restrand_invert=TRUE)
#' nmatlist2heatmaps(nmatlist2i,
#'    title="Test data, centered on inverted tss signal",
#'    transform=c("log2signed", "sqrt"));
#'
#' nmatlist2 <- recenter_nmatlist(nmatlist, recenter_heatmap=1, recenter_range=1000)
#' nmatlist2heatmaps(nmatlist2,
#'    title="Test data, centered on tss signal",
#'    transform=c("log2signed", "sqrt"));
#'
#' nmatlist4 <- restrand_nmatlist(nmatlist2, restrand_invert=FALSE,
#'    transform=c("log2signed"), restrand_heatmap=1,
#'    restrand_range=600, restrand_buffer=100)
#' nmhm <- nmatlist2heatmaps(nmatlist4,
#'    title="Test data, centered on tss signal\nRe-stranded on tss signal",
#'    transform=c("log2signed", "sqrt"));
#'    #k_clusters=8, k_method="euclidean", min_rows_per_k=10, k_heatmap=1,
#'
#' @export
restrand_nmatlist <- function
(nmatlist,
 restrand_heatmap=1,
 restrand_range=NULL,
 restrand_buffer=NULL,
 restrand_invert=FALSE,
 empty_value=0,
 transform=NULL,
 verbose=FALSE,
 ...)
{
   ## Basic steps:

   ## Define the nmatlist with appropriate range to use for recentering
   if (length(restrand_invert) == 0) {
      restrand_invert <- FALSE;
   }
   restrand_invert <- rep(restrand_invert,
      length.out=length(restrand_heatmap));
   if (verbose) {
      if (length(restrand_range) == 1 && restrand_range > 0) {
         jamba::printDebug("restrand_nmatlist(): ",
            "restrand_range: ", restrand_range);
      }
      if (length(restrand_buffer) == 1 && restrand_buffer > 0) {
         jamba::printDebug("restrand_nmatlist(): ",
            "restrand_buffer: ", restrand_buffer);
      }
   }

   # handle optional transformation
   if (length(transform) > 0) {
      transform <- rep(transform, length.out=length(restrand_heatmap));
      transform <- get_numeric_transform(transform)
      if (!is.list(transform)) {
         transform <- list(transform)
      }
   }

   restrand_nmatlist <- lapply(nmatlist[restrand_heatmap], function(nmat){
      if (length(restrand_range) == 1 && restrand_range > 0) {
         nmat <- zoom_nmat(nmat,
            upstream_length=restrand_range,
            downstream_length=restrand_range,
            ...)
      }
      if (length(restrand_buffer) == 1 && restrand_buffer > 0) {
         nmat_buffer <- zoom_nmat(nmat,
            upstream_length=restrand_buffer,
            downstream_length=restrand_buffer,
            ...)
         nmat[, colnames(nmat_buffer)] <- 0;
      }
      nmat;
   })
   nmatlist_nrow <- sapply(nmatlist, nrow);
   if (length(unique(nmatlist_nrow)) > 1) {
      stop("nrow must be identical across all nmatlist entries.")
   }

   ## Define which colnames to use
   restrand_ncols <- lapply(restrand_nmatlist, ncol)
   restrand_colnames_list <- lapply(restrand_nmatlist, colnames)
   # use the superset of colnames available
   restrand_colnames <- jamba::mixedSort(
      Reduce("union", restrand_colnames_list),
      keepNegative=TRUE);
   # expand colnames as needed
   if (length(unique(restrand_ncols)) > 1) {
      jamba::printDebug("restrand_nmatlist(): ",
         "Note that nmatlist contains inconsitent columns.",
         " Adjusting accordingly.");
   }
   # Actually, convert all to numeric matrix for convenience
   restrand_nmatlist <- lapply(seq_along(restrand_nmatlist), function(inum){
      nmat <- restrand_nmatlist[[inum]];
      match_col <- match(restrand_colnames, colnames(nmat));
      nmat_rownames <- rownames(nmat);
      nmat <- jamba::rmNA(naValue=0,
         nmat[, match_col, drop=FALSE]);
      colnames(nmat) <- restrand_colnames;
      rownames(nmat) <- nmat_rownames;
      if (length(transform) > 0) {
         nmat <- transform[[inum]](nmat);
      }
      if (TRUE %in% restrand_invert[[inum]]) {
         nmat <- max(nmat, na.rm=TRUE) - nmat;
      }
      nmat;
   });
   restrand_mat <- Reduce("+", restrand_nmatlist);

   ## Phase One: Determine side with higher coverage per row
   use_rows <- rownames(restrand_mat);
   len1 <- floor(ncol(restrand_mat) / 2);
   cols1 <- seq(from=1, length.out=len1);
   cols2 <- seq(to=ncol(restrand_mat), length.out=len1);
   sum1 <- rowSums(restrand_mat[, cols1], na.rm=TRUE)
   sum2 <- rowSums(restrand_mat[, cols2], na.rm=TRUE)
   use_left <- (sum1 > sum2);
   if (verbose > 1) {
      jamba::printDebug("restrand_nmatlist(): ",
         "data.frame(sum1, sum2, use_left):");
      print(data.frame(sum1, sum2, use_left));
   }

   ## Phase Two: Flip coverage where relevant
   new_nmatlist <- lapply(nmatlist, function(nmat){
      keep12rev <- rev(colnames(nmat))
      new_nmat <- nmat;
      new_nmat[use_left, ] <- nmat[use_left, keep12rev, drop=FALSE]
      attributes(new_nmat) <- attributes(nmat)
      attr(new_nmat, "restrand") <- use_left;
      new_nmat;
   })
   new_nmatlist
}

#' Quick summary of normalizedMatrix data
#'
#' Quick summary of normalizedMatrix data
#'
#' @param nmat `normalizedMatrix` object
#' @param ... additional arguments are ignored.
#'
#' @family jam utility functions
#'
#' @returns `data.frame` with columns:
#'    * `"row"` with rownames,
#'    * `"summit_name"` when `recenter_nmatlist()` has been applied,
#'    * `"restrand"` when `restrand_nmatlist()` has been applied.
#'
#' @export
nmat_summary <- function
(nmat,
 ...)
{
   nmat_df <- data.frame(check.names=FALSE,
      row=attr(nmat, "dimnames")[[1]]);
   if ("summit_name" %in% names(attributes(nmat))) {
      nmat_df$summit_name <- attr(nmat, "summit_name")
   }
   if ("restrand" %in% names(attributes(nmat))) {
      nmat_df$restrand <- attr(nmat, "restrand")
   }
   return(nmat_df);
}


#' Quick summary of a list of normalizedMatrix data
#'
#' Quick summary of a list of normalizedMatrix data
#'
#' @param nmatlist `list` of `normalizedMatrix` objects
#' @param ... additional arguments are ignored.
#'
#' @family jam utility functions
#'
#' @returns `list` of `data.frame` with columns:
#'    * `"row"` with rownames,
#'    * `"summit_name"` when `recenter_nmatlist()` has been applied,
#'    * `"restrand"` when `restrand_nmatlist()` has been applied.
#'
#' @export
nmatlist_summary <- function
(nmatlist,
 ...)
{
   lapply(nmatlist, nmat_summary);
}
