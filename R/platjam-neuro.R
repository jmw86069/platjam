
#' Convert frequency matrix to normalizedMatrix format
#'
#' Convert frequency matrix to normalizedMatrix format
#'
#' This function takes input data with frequency represented
#' as columns, and observations as rows. It will optionally
#' scale each row to have fixed minimum-to-maximum value
#' range, given a range of frequencies to use for scaling.
#'
#' @return `normalizedMatrix` object, a subclass of `matrix`,
#'    as defined in `EnrichedHeatmap`.
#'
#' @param mat numeric matrix with frequency represented as
#'    columns, and whose numeric frequency is stored in
#'    `colnames(mat)` as character values.
#' @param target_frequency numeric vector representing a
#'    range of frequencies to be considered the "target",
#'    and thus highlighted by `EnrichedHeatmap::EnrichedHeatmap()`.
#'    By default the target range is used to order rows from
#'    high to low signal, unless the order is specified
#'    otherwise. When `target_frequency` is `NULL`, there
#'    is no target indicated in the normalizedMatrix.
#' @param ... additional arguments are ignored.
#'
#' @family jam import functions
#'
#' @export
frequency_matrix2nmat <- function
(mat,
 target_frequency=NULL,
 target_name="target",
 signal_name="Frequency",
 do_scale=TRUE,
 scale_from=0,
 scale_to=1,
 scale_frequency=NULL,
 apply_floor=TRUE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to convert a frequency matrix to normalizedMatrix
   x_freq <- as.numeric(colnames(x));
   if (any(is.na(x_freq))) {
      stop("colnames(x) must represent numeric frequency values.");
   }
   #
   freq_index <- seq_len(ncol(x));
   target_index <- which(x_freq >= min(target_frequency) &
         x_freq <= max(target_frequency));
   if (length(target_index) > 0) {
      upstream_index <- which(freq_index < min(target_index));
   } else {
      upstream_index <- NULL;
   }
   downstream_index <- which(freq_index > max(c(0, target_index)));
   if (verbose) {
      printDebug("frequency_matrix2nmat(): ",
         "upstream_index:",
         upstream_index,
         ", target_index:",
         target_index,
         ", downstream_index:",
         downstream_index);
   }

   ## Optionally run normScale() on each row
   if (do_scale) {
      if (length(scale_frequency) > 0) {
         scale_index <- which(x_freq >= min(scale_frequency) &
               x_freq <= max(scale_frequency));
      } else {
         scale_index <- freq_index;
      }
      if (verbose) {
         printDebug("frequency_matrix2nmat(): ",
            "scale_index:",
            scale_index);
      }
      mat <- rowNormScale(mat,
         from=scale_from,
         to=scale_to,
         col_range=scale_index);
      attr(mat, "scale_index") <- scale_index;
      attr(mat, "scale_range") <- c(scale_from,
         scale_to);
      if (apply_floor) {
         if (verbose) {
            scale_diff <- (scale_to - scale_from) / 2;
            printDebug("frequency_matrix2nmat(): ",
               "scale_diff:",
               scale_diff);
         }
         attr(mat, "scale_floor") <- c(scale_from - scale_diff,
            scale_to + scale_diff);
         mat <- noiseFloor(mat,
            minimum=scale_from - scale_diff,
            ceiling=scale_to + scale_diff);
      }
   }

   attr(mat, "upstream_index") <- upstream_index;
   attr(mat, "target_index") <- target_index;
   attr(mat, "downstream_index") <- downstream_index;
   attr(mat, "extend") <- c(0, max(x_freq));
   attr(mat, "smooth") <- FALSE;
   attr(mat, "signal_name") <- signal_name;
   attr(mat, "target_name") <- target_name;
   attr(mat, "target_is_single_point") <- FALSE;
   attr(mat, "background") <- 0;
   attr(mat, "signal_is_categorical") <- FALSE;
   class(mat) <- c("normalizedMatrix", "matrix");
   return(mat);
}


#' Normalize and scale data per row
#'
#' Normalize and scale data per row
#'
#' This function essentially calls `jamba::normScale()`
#' on each row of a numeric matrix. By default, it scales
#' values to a fixed numeric range from 0 to 1, where
#' the minimum value is set to 0 and the maximum value is
#' set to 1. It is much more configurable, see `jamba::normScale()`
#' help docs.
#'
#' @param x,from,to,naValue,singletMethod arguments passed to
#'    `jamba::normScale()`. Note that the default `low,high` values
#'    use the column range defined by `col_range`.
#' @param low numeric value or `NULL`, passed to `jamba::normScale()`
#'    for each row in `x`. When `low` is `NULL`, it uses the
#'    minimum value in `col_range` for each row.
#' @param high numeric value or `NULL`, passed to `jamba::normScale()`
#'    for each row in `x`. When `high` is `NULL`, it uses the
#'    maximum value in `col_range` for each row.
#' @param col_range integer vector referring to column numbers in
#'    the input `x` matrix to use in `jamba::normScale()`. When
#'    `col_range` is `NULL`, it uses all columns in `x`.
#' @param ... additional arguments are passed to `jamba::normScale()`.
#'
#' @family jam utility functions
#'
#' @examples
#' m <- matrix(1:9, ncol=3);
#' m;
#'
#' rowNormScale(m);
#' rowNormScale(m, from=0, to=10);
#'
#' @export
rowNormScale <- function
(x,
 from=0,
 to=1,
 naValue=NA,
 low=NULL,
 high=NULL,
 singletMethod="min",
 col_range=NULL,
 ...)
{
   if (length(col_range) == 0) {
      col_range <- seq_len(ncol(x));
   }
   t(apply(x, 1, function(i){
      if (length(low) == 0) {
         low <- min(i[col_range], na.rm=TRUE);
      }
      if (length(high) == 0) {
         high <- max(i[col_range], na.rm=TRUE);
      }
      jamba::normScale(i,
         low=low,
         high=high,
         from=from,
         to=to,
         singletMethod=singletMethod,
         naValue=naValue,
         ...);
   }));
}

