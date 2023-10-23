#
# platjam-coverhm.R
#
# functions to make it easy to create sophisticated coverage heatmaps
#
# Ultimately a helper function will call nmatlist2heatmaps()
# to make the customization easier for non-R-programmers.


#' Make coverage heatmaps
#'
#' (IN DEV) Make coverage heatmaps using a simple set of config options
#'
#' This function is still in development and is not yet active.
#'
#' This function is intended as a wrapper function which calls
#' `nmatlist2heatmaps()` in a more organized way.
#'
#' The `config_df` is a `data.frame` with one row per coverage heatmap
#' to produce, and all options associated with that heatmap
#' are stored on the same row of the `data.frame`.
#'
#' * file: `character` file path to the coverage matrix file
#' * name: `character` string used as a name, used in difference calculations
#' * label: `character` string used as a label above each heatmap
#' * display: `logical` indicating whether to display each heatmap
#' * group: `character` string indicating a panel group, used to share
#' certain settings across groups of heatmap panels. When this value
#' is not defined, each heatmap is defined as its own `group`.
#' * color: `character` string that refers to a specific color gradient,
#' or a comma-delimited series of colors to use as a color gradient.
#' See below for details.
#' * ceiling: `numeric` used to define the maximum numeric value
#' applied to the color gradient, default=0.7
#' * ylim: `numeric` used to define a fixed y-axis range used for the
#' profile plot atop each heatmap. When this value is absent or `NA`
#' the maximum y-axis value for the `group` is used.
#' * control: `character` string that defines an optional control `name`
#' used to subtract coverage from this row.
#'
#' ## Colors
#'
#' Colors can be define one of a few ways:
#'
#' - name of a color gradient
#' - name of a single color
#' - comma-delimited colors, either as color names, or as hex colors
#' in the format `#FFAABB`.
#'
#' Color gradients from RColorBrewer are recognized. You can view
#' * linear color gradients: `RColorBrewer::display.brewer.all(type="seq")`
#' * divergent color gradients: `RColorBrewer::display.brewer.all(type="div")`
#'
#' Color gradients are also recognized from viridis, however these are
#' only linear: `viridis::inferno()`, `viridis::magma()`, `viridis::plasma()`
#' `viridis::cividis()`,`viridis::viridis()`.
#'
#' @family jam coverage heatmap functions
#'
#' @examples
#' # display RColorBrewer linear gradients
#' RColorBrewer::display.brewer.all(type="seq")
#'
#' # display RColorBrewer divergent gradients
#' RColorBrewer::display.brewer.all(type="div")
#'
#' # viridis linear gradients
#' jamba::showColors(list(
#'    inferno=viridis::inferno(11),
#'    magma=viridis::magma(11),
#'    plasma=viridis::plasma(11),
#'    cividis=viridis::cividis(11),
#'    viridis=viridis::viridis(11)))
#'
#' # jam_linear
#' jamba::showColors(jam_linear)
#'
#' # jam_divergent
#' jamba::showColors(jam_divergent)
#'
# make_coverage_heatmaps <- function
# (config_df,
#  anno_df=NULL,
#  ...)
# {
#    # TODO;
# }



#' Row order for nmatlist coverage heatmap
#'
#' Row order for nmatlist coverage heatmap
#'
#' This function is a simple wrapper function to return the
#' row order for the output of `nmatlist2heatmaps()`.
#' It traverses the `ComplexHeatmap::Heatmap` object,
#' including optional row slices when the rows are
#' partitioned. It returns the actual character
#' vector of rownames, optionally split into a `list`.
#'
#' When the heatmap is fully drawn, for example with
#' `nmatlist2heatmaps(..., do_plot=TRUE)`, the output
#' object includes an element `"HM_drawn"` which also contains
#' the row order as displayed in that heatmap. Otherwise,
#' the output contains element `"draw"` with the heatmap
#' that would be drawn by `ComplexHeatmap::draw()`. In
#' this case the row order has not yet been defined,
#' however this function will evaluate the relevant
#' criteria to determine the row order that would be
#' rendered. And note that this process takes slightly
#' more time, but much less time than rendering the
#' entire set of heatmaps.
#'
#' @family jam coverage heatmap functions
#'
#' @return `list` of rownames, where each list represents one
#'    `row_slice` used in the heatmap, as defined by
#'    `nmatlist2heatmaps()` arguments `partition` or
#'    `k_clusters`.
#'
#' @param HM `list` output from `nmatlist2heatmaps()` that
#'    is expected to contain at least one list element
#'    named `"HM_drawn"` or `"draw"`, with `"HM_drawn"`
#'    used first if present. Either object is also expected
#'    to contain slotName ht_list which represents the
#'    HeatmapList. Among the heatmaps in the HeatmapList,
#'    the first with class `"Heatmap"` that is not named
#'    "cluster" will be used, since "cluster" is the name
#'    of the heatmap used to represent k-means or other
#'    clustering output to partition rows.
#' @param ... additional arguments are ignored.
#'
#' @export
nmathm_row_order <- function
(HM,
 ...)
{
   # HM <- nmathm;
   if ("HM_drawn" %in% names(HM)) {
      HM <- HM[["HM_drawn"]];
   } else if ("draw" %in% names(HM)) {
      HM <- HM$draw$HM_temp;
   }

   if (!"ht_list" %in% slotNames(HM)) {
      stop("Input HM does not contain 'ht_list'.");
   }

   hmsdim <- sdim(HM@ht_list);
   hmwhich <- which(hmsdim$class %in% "Heatmap" & !rownames(hmsdim) %in% "cluster");
   hmwhich1 <- head(hmwhich, 1);
   row_order <- multienrichjam::heatmap_row_order(HM@ht_list[[hmwhich1]]);
   return(row_order);
}


#' Zoom the x-axis range for a list of normalizedMatrix coverage data
#'
#' Zoom the x-axis range for a list of normalizedMatrix coverage data
#'
#' This function filters the matrix columns by distance, and updates
#' important associated attributes:
#'
#' * `attr(nmat, "upstream_index")` - the column index positions upstream the target region
#' * `attr(nmat, "downstream_index")` - the column index positions downstream the target region
#' * `attr(nmat, "target_index")` - the column index positions representing the target region
#' * `attr(nmat, "extend")` - the genomic distance upstream and downstream the target region
#'
#' @family jam coverage heatmap functions
#'
#' @param nmatlist `list` of `normalizedMatrix` objects. Each
#'    `normalizedMatrix` is passed to `zoom_nmat()`.
#' @param upstream_length,downstream_length `numeric` vector whose
#'    values are recycled to length `length(nmatlist)`. Each value is
#'    passed to `zoom_nmat()` so each matrix can be zoomed to independent
#'    ranges.
#' @param ... additional arguments are passed to `zoom_nmat()`.
#'
#' @export
zoom_nmatlist <- function
(nmatlist,
 upstream_length=500,
 downstream_length=500,
 ...)
{
   #
   upstream_length <- rep(
      jamba::rmNULL(upstream_length, nullValue=NA),
      length.out=length(nmatlist));
   downstream_length <- rep(
      jamba::rmNULL(downstream_length, nullValue=NA),
      length.out=length(nmatlist));
   new_nmatlist <- lapply(seq_along(nmatlist), function(inmat){
      zoom_nmat(nmatlist[[inmat]],
         upstream_length=upstream_length[[inmat]],
         downstream_length=downstream_length[[inmat]],
         ...)
   })
   names(new_nmatlist) <- names(nmatlist);
   return(new_nmatlist);
}


#' Zoom the x-axis range for a normalizedMatrix coverage data
#'
#' Zoom the x-axis range for a normalizedMatrix coverage data
#'
#' This function is typically called by `zoom_nmatlist()` but can
#' be called on an individual `normalizedMatrix` object.
#'
#' This function filters the matrix columns by distance, and updates
#' important associated attributes:
#'
#' * `attr(nmat, "upstream_index")` - the column index positions upstream the target region
#' * `attr(nmat, "downstream_index")` - the column index positions downstream the target region
#' * `attr(nmat, "target_index")` - the column index positions representing the target region
#' * `attr(nmat, "extend")` - the genomic distance upstream and downstream the target region
#'
#' @family jam coverage heatmap functions
#'
#' @param nmat `normalizedMatrix` object, where the length extended from
#'    the target region is stored in `attr(nmat, "extend")` as a two-element
#'    integer vector representing upstream, and downstream length.
#'    Each column indicated in `attr(nmat, "upstream_index")` is expected
#'    to represent equal-sized bins spanning that range. Columns are
#'    retained if the farthest distance of the column is less
#'    than `upstream_length`.
#' @param upstream_length,downstream_length `numeric` coordinate maximum
#'    range from the target center region. When either value is `NULL`
#'    no threshold is applied, which is equivalent to `Inf`.
#'    The values are forced positive `abs(upstream_length)` as these
#'    are absolute magnitude length from the target region.
#' @param ... additional arguments are ignored.
#'
#' @export
zoom_nmat <- function
(nmat,
 upstream_length=500,
 downstream_length=500,
 ...)
{
   #
   if (length(upstream_length) == 0 || all(is.na(upstream_length))) {
      upstream_length <- Inf;
   }
   upstream_length <- abs(upstream_length);
   if (length(downstream_length) == 0 || all(is.na(downstream_length))) {
      downstream_length <- Inf;
   }
   downstream_length <- abs(downstream_length);

   # detect bin size
   bin_width <- NULL;
   if (length(attr(nmat, "upstream_index")) > 0) {
      bin_width <- round(attr(nmat, "extend")[1] /  length(attr(nmat, "upstream_index")))
      upstream_start <- rev(seq_along(attr(nmat, "upstream_index")) * bin_width);
   } else {
      upstream_start <- NULL;
   }
   upstream_keep <- (upstream_start <= upstream_length)
   new_extend1 <- max(upstream_start[upstream_keep], na.rm=TRUE);

   if (length(attr(nmat, "downstream_index")) > 0) {
      if (length(bin_width) == 0) {
         bin_width <- round(attr(nmat, "extend")[2] /  length(attr(nmat, "downstream_index")))
      }
      downstream_end <- seq_along(attr(nmat, "downstream_index")) * bin_width;
   } else {
      downstream_end <- NULL;
   }
   downstream_keep <- (downstream_end <= downstream_length)
   new_extend2 <- max(downstream_end[downstream_keep], na.rm=TRUE);
   new_extend <- c(new_extend1, new_extend2);

   target_keep <- rep(TRUE, length(attr(nmat, "target_index")));

   column_set <- c(attr(nmat, "upstream_index"),
      attr(nmat, "target_index"),
      attr(nmat, "downstream_index"))
   column_keep <- column_set[c(upstream_keep,
      target_keep,
      downstream_keep)];

   new_nmat <- nmat[, column_keep, drop=FALSE];
   attr_keep <- setdiff(names(attributes(nmat)),
      c("dim", "dimnames"));
   attributes(new_nmat)[attr_keep] <- attributes(nmat)[attr_keep];
   attr(new_nmat, "upstream_index") <- seq_len(sum(upstream_keep));
   attr(new_nmat, "downstream_index") <- tail(seq_len(ncol(new_nmat)), sum(downstream_keep));
   attr(new_nmat, "target_index") <- setdiff(seq_len(ncol(new_nmat)),
      c(attr(new_nmat, "upstream_index"),
         attr(new_nmat, "downstream_index")));
   attr(new_nmat, "extend") <- new_extend;
   return(new_nmat);
}

