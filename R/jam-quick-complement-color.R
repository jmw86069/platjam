
#' Quick color complement function
#'
#' @family jam color functions
#'
#' @returns `character` vector of colors
#'
#' @examples
#' x <- "dodgerblue";
#' jamba::showColors(list(x=x,
#'    `quick_complement_color(x)`=quick_complement_color(x),
#'    `colorjam::color_complement(x, preset="ryb")`=colorjam::color_complement(x, preset="ryb")))
#'
#' x <- rainbow_hcl(10);
#' jamba::showColors(list(x=x,
#'    `quick_complement_color(x)`=quick_complement_color(x),
#'    `colorjam::color_complement(x, preset="ryb")`=colorjam::color_complement(x, preset="ryb")))
#'
#' @export
quick_complement_color <- function
(x,
 space=c("hls", "hcl"),
 preset="ryb",
 do_plot=FALSE,
 ...)
{
   #
   space <- match.arg(space);

   flip_hue <- function(h, preset="ryb", ...) {
      # x_hw <- colorjam::hw2h(h, ...);
      # x_hw <- colorjam::h2hw(h, ...);
      x_hw <- h;
      x_hw_180 <- (x_hw + 180) %% 360;
      # x_h_180 <- colorjam::h2hw(x_hw_180, ...);
      x_h_180 <- colorjam::hw2h(x_hw_180,
         preset=preset,
         ...);
      return(x_h_180);
   }
   ret_list <- list(input=x);
   # apply using HCL color space
   if ("hcl" %in% space) {
      x_hcl <- jamba::col2hcl(x);
      x_h <- x_hcl["H",];
      x_h_180 <- flip_hue(x_hcl["H",],
         preset=preset,
         ...)
      x_180 <- jamba::hcl2col(H=x_h_180,
         C=x_hcl["C",],
         L=x_hcl["L",]);
      ret_list$output_hcl <- x_180;
   }

   # apply using HLS color space
   if ("hls" %in% space) {
      x_hex <- jamba::rgb2col(col2rgb(x));
      x_rgb <- colorspace::hex2RGB(x_hex);
      x_hls <- t(as(x_rgb, "HLS")@coords);
      x_h_180_2 <- flip_hue(x_hls["H",],
         preset=preset,
         ...)
      x_hls_180 <- colorspace::hex(
         colorspace::HLS(H=x_h_180_2,
            L=x_hls["L",],
            S=x_hls["S",]));
      ret_list$output_hls <- x_hls_180;
   }
   if (do_plot) {
      jamba::showColors(ret_list)
   }
   return(ret_list[[2]]);
}

