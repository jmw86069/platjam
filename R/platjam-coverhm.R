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
#' applied to the color gradient
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
#' showColors(jam_linear)
#'
#' # jam_divergent
#' showColors(jam_divergent)
#'
#' @export
make_coverage_heatmaps <- function
(config_df,
 anno_df=NULL,
 ...)
{
}



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
