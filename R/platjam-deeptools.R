#' Import deepTools coverage matrix to normalizedMatrix
#'
#' Import deepTools coverage matrix to normalizedMatrix
#'
#' This function imports deepTools matrix data to produce
#' a `normalizedMatrix` object as defined in EnrichedHeatmap.
#'
#' The initial use case is to import binned data across a fixed range,
#' relative to some central feature "target", which is defined
#' using upstream distance, downstream distance, and bins defined
#' with fixed width.
#'
#' Minor modifications were made previously to accommodate scaled
#' coverage bins. More development will be done to handle the
#' column headers properly in future.
#'
#' ## Supported features
#'
#' * Recognizes upstream, downstream, bin size
#' * Recognizes sample labels
#'
#'    * multiple samples are split into a `list` of matrix objects
#'
#' * Recognizes row groups
#'
#'    * Row groups are stored in `attr(x, "anno_df")` of the returned object
#'
#' * Duplicated rownames are adjusted using `jamba::makeNames()`
#'
#'
#' @family jam import functions
#'
#'
#' @export
deepTools_matrix2nmat <- function
(x=NULL,
 filename=NULL,
 signal_name=NULL,
 target_name=NULL,
 background=0,
 smooth=FALSE,
 target_is_single_point=FALSE,
 signal_is_categorical=FALSE,
 mat_grep="[-0-9]+:[-0-9]+",
 upstream_grep="^[-]",
 downstream_grep="^[^-]",
 target_grep="^0$",
 verbose=FALSE,
 ...)
{
   # import YAML using the first line only
   x1 <- readLines(filename, n=1)
   xyaml <- yaml::read_yaml(text=gsub("^@", "", x1))
   samples <- xyaml[["sample_labels"]];
   bodysize <- xyaml[["body"]]; # optional for gene body width?
   binsize <- xyaml[["bin size"]];
   upstream <- xyaml[["upstream"]];
   downstream <- xyaml[["downstream"]];

   sort_regions <- xyaml[["sort regions"]]; # descend
   sort_using <- xyaml[["sort using"]]; # mean

   if (length(target_name) == 0) {
      target_name <- xyaml[["ref point"]];
   }
   signal_name <- samples;

   # row group boundaries (row number at each group change)
   group_boundaries <- xyaml[["group_boundaries"]];
   group_labels <- xyaml[["group_labels"]];

   # clean up group_labels
   group_labels <- gsub("[.](bed)(|[.]gz)$", "", ignore.case=TRUE, group_labels);

   # sample boundaries - column number for each column grouping?
   sample_boundaries <- xyaml[["sample_boundaries"]];

   if (TRUE %in% verbose) {
      jamba::printDebug("deepTools_matrix2nmat(): ",
         "\nsamples:", samples,
         "\nrow groups:", group_labels,
         "\nbinsize:", binsize,
         "\nupstream:", upstream,
         "\nbody:", bodysize,
         "\ndownstream:", downstream)
   }


   x <- data.table::fread(text=readLines(filename),
      skip=1,
      header=FALSE,
      sep="\t")
   rownames(x) <- jamba::makeNames(x[[4]]);
   if (TRUE %in% verbose) {
      jamba::printDebug("deepTools_matrix2nmat(): ",
         "Loaded data.")
   }

   mats <- as.matrix(x[, -1:-6, drop=FALSE])

   mat_n <- ceiling(ncol(mats)/length(samples));

   ## required attributes
   starts <- seq(from=-upstream,
      by=binsize,
      length.out=mat_n);
   ends <- starts + binsize - 1;
   upstream_index <- which(starts < 0);
   if (any(bodysize > 0)) {
      target_index <- which(starts == 0);
      downstream_index <- which(starts > 0);
      new_colnames <- paste0("u", seq_len(ncol(im)));
   } else {
      target_index <- numeric(0);
      downstream_index <- which(starts >= 0);
      new_colnames <- paste0(starts, ":", ends);
   }
   # target_name <- xyaml$group_labels;

   # column groupings (optional)
   mat_split <- rep(
      rep(samples, each=mat_n),
      length.out=ncol(mats))
   mat_split <- factor(mat_split, levels=unique(samples));

   # make list of columns in group order as they appear
   colnames_l <- split(colnames(mats), mat_split)

   # split columns to individual matrices
   mat_l <- lapply(seq_along(colnames_l), function(k) {
      i <- colnames_l[[k]];
      im <- mats[, i, drop=FALSE]
      colnames(im) <- new_colnames;
      rownames(im) <- rownames(x)
      attr(im, "signal_name") <- samples[k]
      attr(im, "target_name") <- target_name
      attr(im, "upstream_index") <- upstream_index
      attr(im, "target_index") <- target_index
      attr(im, "downstream_index") <- downstream_index
      attr(im, "extend") <- c(upstream, downstream)
      attr(im, "target_is_single_point") <- target_is_single_point;
      attr(im, "signal_is_categorical") <- signal_is_categorical;

      im
   })

   # add anno_df if relevant
   if (length(group_boundaries) > 2) {
      group_sizes <- diff(group_boundaries)
      anno_df <- data.frame(check.names=FALSE,
         row.names=rownames(mat_l[[1]]),
         # name=rownames(mat_l[[1]]),
         group=rep(group_labels, group_sizes));
      attr(mat_l, "anno_df") <- anno_df;
   }
   attr(mat_l, "rownames") <- rownames(mat_l[[1]]);
   attr(mat_l, "colnames") <- new_colnames;

   mat_l;
}
