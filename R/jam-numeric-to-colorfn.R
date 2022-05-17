
#' Assign color function to a numeric vector
#'
#' Assign color function to a numeric vector
#'
#' This function is called internally by `design2colors()`.
#'
#' @param x `numeric` or `integer` vector
#' @param color `character` color
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
#' @param ... additional arguments are ignored.
#'
#' @export
assign_numeric_colors <- function
(x,
 color,
 ratio_for_zero_baseline=3,
 lower_range_expansion=0.2,
 base_color="#fff6f4",
 restrict_pretty_range=TRUE,
 lens=1,
 ...)
{
   #
   irange <- range(x, na.rm=TRUE);
   ispan <- diff(irange, na.rm=TRUE);
   all_integers <- ("integer" %in% class(x) |
         all((x %% 1) == 0));

   if (min(irange) >= 0) {
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
      color1 <- color;
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
#' @param x `list` of `character` color vectors, or `function` defined by
#'    `circlize::colorRamp2()` which encodes rgb colors and breaks
#'    in attributes `"colors"` and `"breaks"` respectively.
#' @param do_print `logical` indicating whether to print output.
#'    The `list` is returned invisibly.
#' @param colorized `logical` indicating whether to print output using
#'    `jamba::printDebug()` with colorized output.
#' @param ... additional arguments are ignored.
#'
#' @export
print_color_list <- function
(x,
 do_print=TRUE,
 colorized=TRUE,
 ...)
{
   x_str <- lapply(x, function(i){
      if (is.function(i)){
         jamba::nameVector(jamba::rgb2col(attr(i, "colors")), attr(i, "breaks"))
      } else {
         i
      }
   });
   if (do_print) {
      if (colorized) {
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
