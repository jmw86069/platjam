
#' Assign color function to a numeric vector
#'
#' Assign color function to a numeric vector
#'
#' This function is called internally by `design2colors()`. It takes
#' a `numeric` vector as input, and applies a color gradient
#' defined by `color` across the range of `numeric` values.
#'
#' When there are negative values, a divergent color scale is generated:
#' * When `color` is a single R color, it is used as the positive color
#' gradient. A complementary color is chosen as the negative color
#' using `quick_complement_color()`. The middle color is `base_color`.
#' * When `color` is a named color gradient, it is assumed to be a
#' divergent gradient. The middle color of the gradient is assigned to zero,
#' and the color gradient is symmetric above and below zero.
#'
#' The `numeric` range is either defined by the data `x`, or when
#' `color_max` is supplied, it is used to define the `numeric`
#' value at which (or higher) the maximum color is assigned.
#'
#' * When `color_max` is not assigned, it is defined as `max(abs(x))`.
#' * If there are negative values, the color gradient will be
#' defined from `-color_max` to `color_max`, with the middle color
#' assigned to zero, as described above.
#' * If there are only positive values, the baseline will be defined
#' using `ratio_for_zero_baseline`, and the maximum color will be defined
#' at `color_max` (and higher).
#'
#' The `ratio_for_zero_baseline` is used only when there are
#' no negative values, since it determines whether the color gradient
#' should be applied starting at zero, or based upon `min(x)`.
#' When the lowest value is less than `ratio_for_zero_baseline * color_max`
#' the color gradient is applied starting from zero;
#' otherwise the baseline uses `min(x)`. In the latter case, non-zero
#' baseline, the lowest value is extended by multiplying
#' `lower_range_expansion` by the span of the numeric range, in order
#' to prevent the lowest non-zero value from being a blank color
#' defined by `base_color`.
#'
#' @family jam color functions
#'
#' @returns `function` as defined by `circlize::colorRamp2()` which takes
#'    a `numeric` vector as input, and returns `character` vector of
#'    colors as output. The attributes described below are used to
#'    show a suitable summary of colors in `jamba::showColors()`;
#'    and are used by `jamses::heatmap_se()` to define a
#'    usable color legend with reasonable number of labels.
#'    The `function` also contains two important attributes:
#'    * `"breaks"`: `numeric` values at break positions used when the
#'    color function was defined. Note the values are taken from the
#'    color function environment, so modifying breaks directly will
#'    not affect the color function output.
#'    * `"colors"`: a matrix of colors in rgb format, with columns
#'    red, blue, green, and `numeric` values ranging from 0 to 1.
#'    The number of rows equals the length of color breaks.
#'
#' @param x `numeric` or `integer` vector
#' @param color `character` color
#' @param color_max `numeric` optional fixed value to define the
#'    `numeric` value associated with the maximum color gradient color.
#' @param ratio_for_zero_baseline `numeric` indicating the ratio of
#'    max to min numeric range, above which the baseline value
#'    should include zero. This argument is intended when values
#'    include `c(3, 100)` and the numeric range is more reasonable
#'    to include zero; compared to values between `c(185, 192)` where
#'    it makes more sense to apply the color gradient only near
#'    those values, and not to use `c(0, 192)`.
#' @param lower_range_expansion `numeric` used when the baseline is not
#'    zero, see `ratio_for_zero_baseline`, to adjust the baseline below
#'    the lowest value by this fraction of the span of values (max - min)
#'    below the lowest value. The intention is to prevent assigning
#'    white `base_color` to this lowest value, instead it should use
#'    a very pale variant of `color`.
#' @param base_color `character` string used as the baseline color when
#'    defining a color gradient to `color`. This value is passed to
#'    `jamba::getColorRamp(..., defaultBaseColor)`.
#' @param restrict_pretty_range `logical` indicating whether the values
#'    returned by `pretty()` should be filtered by `range(x)`, excluding
#'    any numeric breaks outside that range. The practical effect is that
#'    it reduces labels displayed by `jamba::showColors()` or by color
#'    legends drawn by `ComplexHeatmap::Heatmap()`.
#' @param lens `numeric` value that defines the intensity of color gradient,
#'    with higher values increasing the rate of increase in color intensity,
#'    and negative values decreasing this rate. This argument is passed
#'    to `jamba::getColorRamp()`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
assign_numeric_colors <- function
(x,
 color,
 color_max=NULL,
 ratio_for_zero_baseline=3,
 lower_range_expansion=0.2,
 base_color="#fff6f4",
 restrict_pretty_range=TRUE,
 lens=1,
 verbose=FALSE,
 ...)
{
   if (inherits(x, "integer64")) {
      x <- as.numeric(x)
   }
   # define numeric range of the data
   irange <- range(x, na.rm=TRUE);

   # optionally apply numeric ceiling
   color_max <- head(jamba::rmNA(color_max), 1);
   if (length(color_max) > 0) {
      irange <- range(
         jamba::noiseFloor(c(irange, color_max),
            minimum=-abs(color_max),
            ceiling=abs(color_max)))
   }

   # define spanned numeric range
   ispan <- diff(irange, na.rm=TRUE);
   all_integers <- ("integer" %in% class(x) |
         all((x %% 1) == 0));

   if (min(irange) >= 0) {
      if (verbose) {
         jamba::printDebug("assign_numeric_colors(): ",
            "Values are positive, defining a linear gradient.")
      }
      # use a linear color scale if the minimum value is at least zero
      if (min(irange) > 0 && (
         # if lowest value is relatively close to zero, use baseline zero
         max(irange) / min(irange) >= ratio_for_zero_baseline ||
            max(irange) == min(irange))) {
         irange[1] <- 0;
      } else if (min(irange) > 0 &&
            length(unique(irange)) > 1) {
         # expand lower end of numeric range to prevent white
         irange <- irange + c(diff(irange) * -lower_range_expansion, 0);
      }

      # linear color scale
      if (length(unique(irange)) == 1) {
         if (irange[1] == 0) {
            irange <- c(0, 1);
            pretty_range <- c(0, 1);
         } else {
            irange <- c(0, irange[2]);
            pretty_range <- c(0, irange[2]);
         }
      } else {
         pretty_range <- pretty(irange);
         # If values are all integers, only use integer pretty_range values.
         # This step fixes when values are c(0, 1), or c(1, 2), the
         # color key only needs to indicate these two values.
         if (all_integers && !all((pretty_range %% 1) == 0)) {
            pretty_range <- pretty_range[(pretty_range %% 1) == 0];
         }
         # restrict pretty_range to observed range
         if (restrict_pretty_range) {
            pretty_range <- pretty_range[pretty_range >= irange[1] &
                  pretty_range <= irange[2]];
         }
         if (min(pretty_range) > 0) {
            # if lowest value is relatively close to zero, use baseline zero
            if (max(pretty_range) / min(pretty_range) >= ratio_for_zero_baseline) {
               pretty_range <- c(0, pretty_range);
            }
         }
      }
      #color1 <- jamba::color2gradient(color, n=3, dex=10);
      # one color function with specific color breaks
      icolors1 <- circlize::colorRamp2(
         breaks=seq(from=irange[1],
            to=irange[2],
            length.out=6),
         colors=jamba::getColorRamp(color,
            defaultBaseColor=base_color,
            #c(base_color, color1),
            lens=lens,
            n=6));
      # use that color function with pretty_breaks
      icolors <- circlize::colorRamp2(
         breaks=pretty_range,
         colors=icolors1(pretty_range));
   } else {
      if (verbose) {
         jamba::printDebug("assign_numeric_colors(): ",
            "Values are negative, defining a divergent gradient.")
      }
      # divergent color scale
      irange <- max(abs(irange)) * c(-1, 1);
      pretty_range <- pretty(irange);
      # If values are all integers, only use integer pretty_range values.
      # This step fixes when values are c(0, 1), or c(1, 2), the
      # color key only needs to indicate these two values.
      if (all_integers && !all((pretty_range %% 1) == 0)) {
         pretty_range <- pretty_range[(pretty_range %% 1) == 0];
      }
      # restrict pretty_range to observed range
      if (restrict_pretty_range) {
         pretty_range <- pretty_range[pretty_range >= irange[1] &
               pretty_range <= irange[2]];
      }

      # define divergent color scale
      color1 <- tryCatch({
         jamba::rgb2col(col2rgb(color))
      }, error=function(e){
         color
      })
      if (verbose) {
         jamba::printDebug("assign_numeric_colors(): ",
            "color1", paste0("'", color1, "'"))
      }

      if (grepl("^#", color1)) {
         color2 <- quick_complement_color(color1);
         # one color function with specific color breaks
         icolors1 <- circlize::colorRamp2(
            breaks=seq(from=irange[1],
               to=irange[2],
               length.out=7),
            colors=jamba::getColorRamp(
               c(color2, base_color, color1),
               divergent=TRUE,
               defaultBaseColor=base_color,
               lens=lens,
               n=7));
      } else {
         # one color function with specific color breaks
         icolors1 <- circlize::colorRamp2(
            breaks=seq(from=irange[1],
               to=irange[2],
               length.out=7),
            colors=jamba::getColorRamp(
               c(color1),
               divergent=TRUE,
               defaultBaseColor=base_color,
               lens=lens,
               n=7));
      }
      # use that color function with pretty_breaks
      icolors <- circlize::colorRamp2(
         breaks=pretty_range,
         colors=icolors1(pretty_range));
   }
   return(icolors);
}


#' Print color list of color vectors or color functions
#'
#' Print color list of color vectors or color functions
#'
#' @family jam color functions
#'
#' @param x `list` of `character` color vectors, or `function` defined by
#'    `circlize::colorRamp2()` which encodes rgb colors and breaks
#'    in attributes `"colors"` and `"breaks"` respectively.
#' @param do_print `logical` indicating whether to print output.
#'    The `list` is returned invisibly.
#' @param colorized `logical` indicating whether to print output using
#'    `jamba::printDebug()` with colorized output.
#' @param digits `numeric` number of digits, passed to `signif()` to limit
#'    the number of digits displayed as a label for color function entries
#'    in `x`.
#' @param max_entries `numeric` maximum number of entries displayed in any
#'    given element of `x`.
#' @param ... additional arguments are ignored.
#'
#' @returns `list` named by `names(x)` of the input data,
#'    each element contains a `character` vector of R colors
#'    named by respective label. The names are derived either directly
#'    from the input color vector, or by the `attr(i, "breaks")`
#'    encoded into the color function by `circlize::colorRamp2()`.
#'
#' @examples
#' # generate example data.frame of annotations
#' n <- 100;
#' set.seed(12);
#' anno_df <- data.frame(
#'    group=sample(c("A", "B", "B"),
#'       size=n,
#'       replace=TRUE),
#'    tss_score=rnorm(n) + 4,
#'    h3k4me1_score=rnorm(n) + 4
#' );
#'
#' # generate a list of colors for the data.frame
#' scl <- design2colors(anno_df,
#'    group_colnames="group",
#'    color_sub=c(tss_score="royalblue",
#'       h3k4me1_score="slateblue"),
#'    plot_type="none")
#'
#' # print the list of colors, color functions, in colorized text
#' print_color_list(scl)
#'
#' @export
print_color_list <- function
(x,
 do_print=TRUE,
 colorized=TRUE,
 digits=4,
 max_entries=Inf,
 ...)
{
   if (length(max_entries) == 0) {
      max_entries <- Inf;
   }
   x_str <- lapply(x, function(i){
      if (is.function(i)){
         icolors <- attr(i, "colors");
         ibreaks <- attr(i, "breaks");
         if (is.numeric(ibreaks)) {
            ibreaks <- signif(ibreaks,
               digits=digits)
         }

         if (any(grepl("#[a-zA-Z0-9]+", icolors))) {
            j <- jamba::nameVector(
               attr(i, "colors"),
               attr(i, "breaks"))
         } else {
            j <- jamba::nameVector(
               jamba::rgb2col(icolors),
               attr(i, "breaks"))
         }
      } else {
         j <- i
      }
      head(j, max_entries)
   })
   if (TRUE %in% do_print) {
      if (TRUE %in% colorized) {
         for (i in names(x_str)) {
            jamba::printDebug(paste0("$", i, ""),
               timeStamp=FALSE,
               comment=FALSE);
            jamba::printDebugI(x_str[[i]],
               timeStamp=FALSE,
               comment=FALSE);
         }
      } else {
         print(x_str);
      }
   }
   return(invisible(x_str));
}


#' Mean hue for a vector
#'
#' Mean hue for a vector
#'
#' @family jam color functions
#'
#' @param x `numeric` angles in degrees representing HCL color hue.
#' @param seed `numeric` used with `set.seed(seed)` for reproducibility.
#' @param ... additional arguments are ignored.
#'
#' @export
mean_hue <- function
(x,
 seed=1,
 ...)
{
   if (length(x) <= 1) {
      return(x)
   }
   if (length(seed) == 1) {
      set.seed(seed);
   }
   x <- x + rnorm(length(x)) * 1e-6;
   xmean <- mean(cos(jamba::deg2rad(x)), na.rm=TRUE);
   ymean <- mean(sin(jamba::deg2rad(x)), na.rm=TRUE);
   degrees <- jamba::rad2deg(atan2(y=ymean, x=xmean)) %% 360;
   radius <- sqrt(xmean^2 + ymean^2);
   attr(degrees, "radius") <- radius;

   return(degrees)
}
