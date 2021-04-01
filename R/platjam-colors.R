#
# platjam-colors.R
#
# Colors specifically designed for use in coverage heatmaps.
# The gradients include linear and divergent color gradients,
# they use multi-hue gradients to expand visible clarity,
# and are color-blind-friendly.


#' Jam Linear Color Gradients
#'
#' Jam Linear Color Gradients that are color-blind-friendly.
#'
#' This data contains linear color gradients generated with
#' a multi-hue technique similar to that used by `RColorBrewer`
#' to expand the visual clarity of a linear single-hue color
#' gradient. Colors were chosed to avoid common color-blindness
#' problems, and to maximize visual differences between
#' color gradients.
#'
#' The colors in `jam_linear` are designed so that colors
#' in `jam_divergent` can be chosen for consistency.
#' For example linear color gradient `"firebrick"` can
#' be used to show coverage values, and `"skyblue_firebrick"`
#' can be used to show coverage difference from control,
#' and in both cases the maximum color is `"firebrick"`.
#'
#' These linear gradients are designed with a white background
#' color. Note the `jam_divergent` colors are designed with
#' a black background, intended to help indicate that these
#' colors are divergent.
#'
#' In general, there are seven warm color tones, and seven cool color
#' tones. Visual distinction is only expected between warm and
#' cool color tones, and is not distinct within the same color
#' tone.
#'
#' Each color is named by the closest corresponding R color:
#'
#' * firebrick
#' * orangered
#' * tomato
#' * sienna
#' * goldenrod
#' * gold
#' * skyblue
#' * dodgerblue
#' * royalblue
#' * slateblue
#' * orchid
#' * magenta
#' * maroon
#'
#' @examples
#' data(jam_linear)
#' jamba::showColors(jam_linear)
#'
#' # show the same with dichromat adjustment
#' jamba::showColors(lapply(jam_linear, dichromat::dichromat))
#'
"jam_linear"


#' Jam Divergent Color Gradients
#'
#' Jam Divergent Color Gradients that are color-blind-friendly.
#'
#' This data contains divergent color gradients generated
#' similar to `jam_linear` with multi-hue technique based upon
#' `RColorBrewer`. Each gradient is designed so that
#' the two colors are visibly distinct even under three
#' different color blindness simulations from `dichromat::dichromat()`.
#'
#' Each gradient is named by the closest corresponding R colors,
#' where the same names are used in `jam_linear` so that colors
#' can be matched where appropriate:
#'
#' * skyblue_firebrick
#' * dodgerblue_orangered
#' * royalblue_tomato
#' * slateblue_sienna
#' * orchid_orange
#' * magenta_goldenrod
#' * maroon_gold
#' * firebrick_skyblue
#' * orangered_dodgerblue
#' * tomato_royalblue
#' * sienna_slateblue
#' * orange_orchid
#' * goldenrod_magenta
#' * gold_maroon
#'
#' These linear gradients are designed with a white background
#' color.
#' These divergent gradients are designed with
#' a black background, intended to help indicate that these
#' colors are divergent. Note that linear gradients `jam_linear`
#' are designed with a white background color.
#'
#' @examples
#' data(jam_divergent)
#' jamba::showColors(jam_divergent)
#'
#' # show the same with dichromat adjustment
#' jamba::showColors(lapply(jam_divergent, dichromat::dichromat))
#'
"jam_divergent"


#' Make divergent color gradient
#'
#' Make divergent color gradient using jam_linear and jam_divergent.
#'
#' @examples
#' jamba::showColors(jam_linear)
#'
#' jg1 <- make_jam_divergent("royalblue", "orangered")
#' jamba::showColors(jg1)
#' showDichromat(jg1)
#'
#' jg2 <- make_jam_divergent("slateblue", "firebrick", n=21)
#' jamba::showColors(jg2)
#' showDichromat(jg2)
#'
#' jg3 <- make_jam_divergent("slateblue", "firebrick", lite=FALSE, n=21)
#' jamba::showColors(jg3)
#' showDichromat(jg3)
#'
#' jg4 <- make_jam_divergent("Blues", "Reds", lite=TRUE, n=21)
#' jamba::showColors(c(jg4,
#'    list(BuRd=jamba::getColorRamp("RdBu_r", n=21))))
#'
#' jg5 <- make_jam_divergent("inferno", "dodgerblue1", lite=FALSE, n=21, gradientWtFactor=1)
#' jamba::showColors(jg5)
#'
#' xseq <- seq(from=-1, to=1, by=0.1);
#' mseq <- matrix(xseq, ncol=1);
#' m <- mseq %*% t(mseq);
#' rownames(m) <- seq_len(nrow(m));
#' colnames(m) <- seq_len(ncol(m));
#' hm1 <- ComplexHeatmap::Heatmap(m[,1:10],
#'    cluster_columns=FALSE,
#'    cluster_rows=FALSE,
#'    row_names_side="left",
#'    border=TRUE,
#'    heatmap_legend_param=list(
#'       border="grey10",
#'       at=seq(from=-1, to=1, by=0.25),
#'       color_bar="discrete"),
#'    col=jg3[[1]])
#'
#' hm2 <- ComplexHeatmap::Heatmap(m[21:1,12:21],
#'    cluster_columns=FALSE,
#'    cluster_rows=FALSE,
#'    border=TRUE,
#'    heatmap_legend_param=list(
#'       border=TRUE,
#'       at=seq(from=-1, to=1, by=0.25),
#'       color_bar="discrete"),
#'    col=jg2[[1]])
#' hm1 + hm2
#'
#' # same as above but showing where to use lens
#' hm3 <- ComplexHeatmap::Heatmap(m[,1:10],
#'    cluster_columns=FALSE,
#'    cluster_rows=FALSE,
#'    row_names_side="left",
#'    border=TRUE,
#'    heatmap_legend_param=list(
#'       border=TRUE,
#'       at=seq(from=-1, to=1, by=0.25),
#'       color_bar="discrete"),
#'    col=jamba::getColorRamp(jg3[[1]], divergent=TRUE, lens=2))
#'
#' hm4 <- ComplexHeatmap::Heatmap(m[21:1,12:21],
#'    cluster_columns=FALSE,
#'    cluster_rows=FALSE,
#'    border=TRUE,
#'    heatmap_legend_param=list(
#'       border=TRUE,
#'       at=seq(from=-1, to=1, by=0.25),
#'       color_bar="discrete"),
#'    col=jamba::getColorRamp(jg2[[1]], divergent=TRUE, lens=2))
#' hm3 + hm4
#'
#' @export
make_jam_divergent <- function
(linear1,
 linear2,
 lite=TRUE,
 n=21,
 ...)
{
   #if (!all(linear1 %in% names(jam_linear))) {
   #   stop("linear1 must be values in names(jam_linear)")
   #}
   #if (!all(linear2 %in% names(jam_linear))) {
   #   stop("linear2 must be values in names(jam_linear)")
   #}
   n_out <- max(c(
      length(linear1),
      length(linear2),
      length(lite),
      length(n)));
   linear1 <- rep(linear1, length.out=n_out);
   linear2 <- rep(linear2, length.out=n_out);
   lite <- rep(lite, length.out=n_out);
   n <- rep(n, length.out=n_out);

   get_jam_gradient <- function(x, lite=TRUE, n=11, ...) {
      if (length(x) == 1 && x %in% names(jam_linear)) {
         if (lite) {
            xcolors <- jam_linear[[x]];
         } else {
            xwhich <- match(x, gsub("_.+", "", names(jam_divergent)));
            xcolors <- jam_divergent[[xwhich]];
            xcolors <- rev(head(xcolors, ceiling(length(xcolors)/2)));
         }
         jamba::getColorRamp(xcolors,
            n=n,
            ...)
      } else {
         if (lite) {
            defaultBaseColor <- "white";
         } else {
            defaultBaseColor <- "black";
         }
         xcolors <- jamba::getColorRamp(x,
            n=n,
            defaultBaseColor=defaultBaseColor,
            ...)
      }
   }

   if (length(names(linear1)) == 0) {
      if (is.atomic(linear1)) {
         names(linear1) <- linear1;
      } else {
         names(linear1) <- jamba::makeNames(rep("linear1", length(linear1)));
      }
   }
   if (length(names(linear2)) == 0) {
      if (is.atomic(linear2)) {
         names(linear2) <- linear2;
      } else {
         names(linear2) <- jamba::makeNames(rep("linear2", length(linear2)));
      }
   }
   gradient_names <- paste0(
      names(linear1),
      "_",
      names(linear2));
   gradient_names <- ifelse(lite,
      gradient_names,
      paste0(gradient_names, "_dark"));
   gradient_list <- lapply(seq_along(gradient_names), function(k){
      nk <- ceiling(n[k] / 2);
      if (nk == 0) {
         nk <- 11;
      }
      gr1 <- head(rev(
         get_jam_gradient(linear1[[k]],
            lite=lite[k],
            n=nk,
            ...)), -1);
      gr2 <- get_jam_gradient(linear2[[k]],
         lite=lite[k],
         n=nk,
         ...);
      if (1 == 2) {
         if (lite[k]) {
            gr1 <- head(rev(jam_linear[[linear1[k]]]), -1)
            gr2 <- jam_linear[[linear2[k]]];
         } else {
            n1 <- jamba::vigrep(paste0("^", linear1[k], "_"), names(jam_divergent));
            n2 <- jamba::vigrep(paste0("_", linear2[k], "$"), names(jam_divergent));
            gr1 <- head(jam_divergent[[n1]], 10)
            gr2 <- tail(jam_divergent[[n2]], 11)
         }
      }
      gr12 <- c(gr1, gr2);
      if (n == 0) {
         gr12 <- jamba::getColorRamp(gr12,
            n=NULL,
            divergent=TRUE)
      }
      gr12;
   })
   names(gradient_list) <- gradient_names;
   return(gradient_list);
}

#' Show colors using dichromat color blindness adjustment
#'
#' Show colors using dichromat color blindness adjustment
#'
#' @param x `list` or `character` vector with R compatible colors.
#' @param type `character` passed to `dichromat::dichromat()` for one
#'    or more types of color blindness to simulate.
#' @param sep `character` used as a delimited to label each resulting
#'    color vector.
#' @param spacer `logical` indicating whether to include a blank spacer
#'    between sets of colors. This spacer is mainly useful for display.
#' @param original `logical` indicating whether to include original colors
#'    and adjusted colors.
#' @param do_plot `logical` indicating whether to plot the results
#'    using `jamba::showColors()`.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' showDichromat(jam_linear["firebrick"])
#'
#' showDichromat(jam_linear[1:2])
#'
#' showDichromat(jam_linear[7:9])
#'
#' showDichromat(jam_linear, type="tritan", spacer=FALSE)
#'
#' showDichromat(jam_linear, type="tritan", spacer=FALSE, original=FALSE)
#'
#' @export
showDichromat <- function
(x,
 type=c("deutan", "protan", "tritan"),
 sep="\n",
 spacer=TRUE,
 original=TRUE,
 do_plot=TRUE,
 ...)
{
   if (is.atomic(x)) {
      x <- list(x);
   }
   if (length(names(x)) == 0) {
      names(x) <- seq_along(x);
   }
   type <- match.arg(type, several.ok=TRUE);
   names(type) <- type;
   x_new <- lapply(seq_along(x), function(i){
      x_di <- lapply(type, function(k){
         dichromat::dichromat(x[[i]], type=k)
      })
      names(x_di) <- paste0(names(x)[i], sep, type);
      if (original) {
         c(x[i], x_di)
      } else {
         x_di
      }
   })
   x_set <- x_new[[1]];
   if (spacer) {
      blank <- list(` `="transparent")
   } else {
      blank <- NULL
   }
   for (i in tail(seq_along(x_new), -1)) {
      x_set <- c(x_set,
         blank,
         #` `=c("transparent"),
         x_new[[i]]);
   }
   #x_set <- unlist(recursive=FALSE, x_new);
   if (do_plot) {
      jamba::showColors(x_set,
         ...)
   }
   return(invisible(x_set));
}

#' Create two-step linear gradient
#'
#' Create two-step linear gradient by gradually blending
#' two linear color gradients
#'
#' This function is intended to produce a two-step linear gradient
#' effect, similar to the strategy used by `RColorBrewer`, but
#' without specific color constraints. See examples.
#'
#' This function takes two color gradients and blends them
#' using a weighting scheme that begins with 100% `color1`, and
#' gradually becomes 100% `color2`.
#'
#' The input `color1` and `color2` can be any input recognized
#' by `jamba::getColorRamp()`. For example a single color can
#' be used to create a gradient, or the name of a known color
#' gradient can be used, for example `"Reds"` will refer
#' to `RColorBrewer` palette `"Reds"`. See the examples.
#'
#' In general most gradients can be blended using this function
#' to produce a new color gradient where both the visual intensity
#' and color hue vary along the gradient, making each color step
#' more visibly distinct than when only the visual intensity
#' changes.
#'
#' When supplying a single color as input to `color1` or `color2`
#' it sometimes works best to alter the brightness of one or both
#' colors so the intermediate gradients have similar intensities.
#' Experimenting with `debug=TRUE` is recommended.
#'
#'
#' @examples
#' ts <- twostep_gradient("yellow", debug=TRUE)
#'
#' ts1 <- twostep_gradient("orange2", "firebrick", n=11, debug=TRUE)
#' ts2 <- twostep_gradient("aquamarine", "dodgerblue", n=11, debug=TRUE)
#'
#' # stitch them together with make_jam_divergent()
#' ts1ts2 <- make_jam_divergent(list(ts2=ts2), list(ts1=ts1), n=21)
#' jamba::showColors(ts1ts2)
#' ts1ts2flat <- make_jam_divergent("dodgerblue", "firebrick", n=21)
#' jamba::showColors(list(
#'    twostep=ts1ts2[[1]],
#'    flat=ts1ts2flat[[1]]))
#'
#' ts3 <- twostep_gradient("Greens", "Blues", n=11, debug=TRUE)
#'
#' ts4 <- twostep_gradient("slateblue2", "firebrick", n=11, debug=TRUE)
#'
#' ts5 <- twostep_gradient("cividis", "inferno", n=11, debug=TRUE, adjust=-1.2)
#'
#' gr1 <- twostep_gradient("slateblue", "purple", debug=TRUE)
#' gr2 <- twostep_gradient("gold", "darkorange", debug=TRUE)
#' div12 <- make_jam_divergent(list(gr1=gr1), list(gr2=gr2))
#' jamba::showColors(div12)
#' div12flat <- make_jam_divergent("purple", "gold")
#' jamba::showColors(list(
#'    twostep=div12[[1]],
#'    flat=div12flat[[1]]))
#'
#' gr1d <- twostep_gradient("slateblue1", "purple", debug=TRUE, lite=FALSE)
#' gr2d <- twostep_gradient("darkorange", "gold", debug=TRUE, lite=FALSE)
#' div12d <- make_jam_divergent(list(gr1d=gr1d), list(gr2d=gr2d))
#' jamba::showColors(div12d)
#' div12dflat <- make_jam_divergent("purple", "gold", lite=FALSE)
#' jamba::showColors(list(
#'    twostep=div12d[[1]],
#'    flat=div12dflat[[1]]))
#'
#' @param color1 `character` color or name of a recognized color gradient.
#' @param color2 `character` color or name of a recognized color gradient;
#'    or when `color2=NULL` then the hue of `color1` is shifted to
#'    emulate the effect of having a similar neighboring color hue.
#'    In this case the input `color1` is used as `color2` to become
#'    the primary output color.
#' @param n `integer` number of gradient colors to return. When `n=0`
#'    or `n=NULL` the output is a color function.
#' @param lite `logical` indicating whether the background color
#'    should be white, or when `lite=FALSE` the background color
#'    is black.
#' @param defaultBaseColor `character` used to define a specific
#'    background color, and therefore overrides `lite`.
#' @param adjust `numeric` value used to adjust the relative
#'    weight between `color1` and `color2`, where values higher
#'    than 1 favor `color2` and negative values, or values less
#'    than 1 favor `color1`.
#' @param do_fixYellow `logical` indicating whether to call
#'    `jamba::fixYellow()` which fixes the greenish hue that
#'    sometimes results from what is intended to be pure yellow.
#' @param debug `logical` indicating whether to create a plot
#'    to show the color blending steps.
#' @param ... additional arguments are passed to `jamba::getColorRamp()`.
#'
#' @export
twostep_gradient <- function
(color1=NULL,
 color2=NULL,
 n=11,
 lite=TRUE,
 defaultBaseColor=NULL,
 adjust=1.5,
 do_fixYellow=TRUE,
 debug=FALSE,
 ...)
{
   nk <- jamba::noiseFloor(n, minimum=5);

   # special case where color2 is NULL
   if (length(color2) == 0) {
      color2 <- color1;
      color1hcl <- jamba::col2hcl(color2);
      H_add <- ifelse(color1hcl["H",] < 40, 50,
         ifelse(color1hcl["H",] > 200, -20,
            ifelse(color1hcl["H",] < 120, -60, 70)));
      color1hcl["H",] <- colorjam::hw2h(preset="ryb",
         colorjam::h2hw(color1hcl["H",], preset="ryb") + H_add);
      color1hcl["L",] <- jamba::noiseFloor(color1hcl["L",],
         minimum=45);
      color1hcl["C",] <- jamba::noiseFloor(color1hcl["C",],
         minimum=160);
      color1 <- jamba::hcl2col(color1hcl)
      showColors(c(color1, color2))
   }
   if (length(defaultBaseColor) == 0) {
      if (lite) {
         defaultBaseColor <- "white"
      } else {
         defaultBaseColor <- "black"
      }
   }
   g1 <- jamba::getColorRamp(color1,
      n=nk,
      defaultBaseColor=defaultBaseColor,
      ...);
   g2 <- jamba::getColorRamp(color2,
      n=nk,
      defaultBaseColor=defaultBaseColor,
      ...);

   # gradient weight
   wseq <- seq(from=1, to=0, length.out=nk - 1);
   if (adjust > 1) {
      wseq <- wseq ^ adjust
   } else if (adjust > 0 && adjust < 1) {
      wseq <- wseq ^ adjust
   } else if (adjust < 0) {
      wseq <- wseq ^ (1/-adjust)
   }
   w1 <- c(1, wseq);
   w2 <- 1 - w1;
   wdf <- data.frame(w1=w1, w2=w2);
   print(wdf);

   g12 <- sapply(seq_len(n), function(i){
      colorjam::blend_colors(c(
         jamba::alpha2col(g1[i], alpha=w1[i]),
         jamba::alpha2col(g2[i], alpha=w2[i])
      ))})

   # optionally "fix" yellow hues
   if (do_fixYellow) {
      g12 <- jamba::fixYellow(g12);
   }

   if (debug) {
      jamba::showColors(list(g1=g1, g2=g2, g12=g12))
      lines(x=seq_len(nk),
         y=-1 * wdf$w1 + 2,
         type="b", lwd=2)
      lines(x=seq_len(nk),
         y=-1 * wdf$w2 + 2,
         type="b", lwd=2)
   }
   if (!n == nk) {
      if (n == 0) {
         n <- NULL
      }
      g12 <- jamba::getColorRamp(g12,
         n=n,
         divergent=TRUE);
   }
   return(g12);
}

