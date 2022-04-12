
#' Convert experiment design into categorical colors
#'
#' Convert experiment design into categorical colors
#'
#' The general goal is to assign categorical colors relevant to
#' the experimental design of an analysis. The basic logic:
#'
#' 1. Assign categorical colors to broadly defined experimental groups.
#' 2. Shade these colors light-to-dark based upon secondary factors.
#' 3. For step 1 above, optionally assign similar color hues by class.
#'
#' ## Categorical colors
#'
#' At its simplest a set of groups can be assigned categorical colors.
#'
#' * colors should be visibly distinct from one another
#' * colors should generally be distinct across forms of color-blindness
#' * colors should be consistent across plots, figures, tables
#'
#' Finally, colors may be pre-defined using a named vector of colors.
#' These colors will be propagated to other entries in the table.
#'
#' ## Light-to-dark gradient
#'
#' The light-to-dark gradient is intended for ordered sub-divisions,
#' for example:
#'
#' * across time points in a time series
#' * across treatment doses in an ordered series
#' * across ordered measurements first-to-last
#'
#' ## Group class
#'
#' The group classification is intended to assign color hues
#' for similar groups:
#'
#' * antagonists, agonists, untreated
#' * treated, untreated
#' * wildtype, mutant form 1, mutant form 2, etc.
#'
#' For example, antagonists may be assigned colors blue-to-purple;
#' agonists may be assigned colors red-to-orange; with a pronounced
#' color hue "gap" between antagonists and agonists.
#'
#' ## Additional categorical color assignment
#'
#' Finally, other annotations associated with samples are assigned
#' categorical colors, visibly distinct from other color assignments.
#'
#' For entries associated with only one design color, for example "Sample_ID",
#' "Sample Name", "Lane Number", or "Well Number",
#' they inherit the design color.
#'
#' For entries associated with more than one design color, for example
#' "Batch", "Date", or perhaps "Dose", they will be assigned a unique
#' color.
#'
#' * additional annotations unique to design colors inherit the design colors
#' * additional categorical colors should not duplicate existing colors
#'
#' ## Future ideas
#'
#' * Assign "additional factors" to colors based upon `class`
#'
#'    * Currently "additional factors" are only tested by class_group and
#'    class_group_lightness.
#'    * It could be useful to test versus `class` alone (if supplied)
#'    * Goal would be to assign color hue using the mean color hue in the class.
#'    * Otherwise the class may be assigned a color inconsistent with the
#'    range of color hues.
#'
#' * Handle numeric columns by applying color gradient
#'
#'    * A truly numeric column (not just integer index values) could
#'    use `circlize::colorRamp2()` to apply color gradient
#'
#'
#' @param x `data.frame` with columns to be colorized
#' @param group_colnames `character` or `intger` vector indicating
#'    which `colnames(x)` to use, in order, for group color assignment.
#' @param lightness_colnames `character` or `intger` vector indicating
#'    which `colnames(x)` to use, in order, for group lightness gradient.
#' @param class_colnames `character` or `intger` vector indicating
#'    higher-level grouping of `group_colnames`
#' @param preset `character` string passed to `colorjam::h2hwOptions()`,
#'    which defines the hues around a color wheel, used when selecting
#'    categorical colors.
#' @param phase,rotate_phase `integer` value, `phase` is passed to
#'    `colorjam::rainbowJam()` to define the light/dark pattern phasing,
#'    which has 6 positions, and negative values reverse the order.
#'    Categorical colors are assigned to the class/group combinations,
#'    after which `phase + rotate_phase` is used for categorical colors
#'    for any remaining values.
#' @param class_pad `integer` zero or greater, indicating the number
#'    of empty hues to insert as a spacer between hues when the class
#'    changes in a sequence of class/group values. Higher values will
#'    ensure the hues in each class are more distinct from each other
#'    across class, and more similar to each other within class.
#' @param end_hue_pad `integer` used to pad hues at the end of a
#'    color wheel sequence, typically useful to ensure the last color
#'    is not similar to the first color.
#' @param desat `numeric` vector extended to length=2, used to desaturate
#'    class/group colors, then remaining colors, in order. The intended
#'    effect is to have class/group colors visibly more colorful than
#'    remaining colors assigned to other factors.
#' @param dex `numeric` vector passed to `jamba::color2gradient()` to
#'    define the darkness expansion factor, where 1 applies a moderate
#'    effect, and higher values apply more dramatic light-to-dark
#'    effect. When `dex` has length=2, the second value is used only
#'    for columns where colors are assigned by `colnames(x)`
#'    using `color_sub`.
#' @param color_sub `character` vector of R colors, where `names(color_sub)`
#'    assign each color to a character string. It is intended to allow
#'    specific color assignments upfront.
#'    * `colnames(x)`: when `names(color_sub)` matches a column name in `x`,
#'    the color is assigned to that color using a color gradient across
#'    the unique character values in that column. Values are assigned in
#'    order of their appearance in `x` unless the column is a `factor`,
#'    in which case colors are assigned to `levels`.
#' @param plot_type `character` string indicating a type of plot for results:
#'    * `"table"`: plots a color table equal to the input `data.frame` where
#'    background cells indicate color assignments.
#'    * `"list"`: plots colors defined for each column in `x` using
#'    `jamba::showColors()`
#'    * `"none"`: no plot is produced
#' @param return_type `character` string indicating the data format to return:
#'    * `"list"`: a `list` of colors, named by `colnames(x)`.
#'    * `"df"`: a `data.frame` in order of `x` with colors assigned to each cell.
#'    * `"vector"`: a `character` vector of R colors, named by assigned
#'    factor level.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to downstream functions.
#'
#' @return output depends upon argument `return_type`:
#' * `"list"`: returns a `list` of colors defined by `colnames(x)`,
#' suitable for use with `ComplexHeatmap::HeatmapAnnotation()` for example.
#' * `"df"`: returns `data.frame` of colors with same dimensions as the
#' input `x`. Suitable for use with `jamba::imageByColors()` for example.
#' * `"vector"`: returns `character` vector of R colors, whose names represent
#' values in `x`, where the values should be substituted with the color.
#' Suitable for use with `ggplot2::color_manual(values=colors)`.
#'
#'    In all cases, the `attributes()` of the returned object also includes
#'    colors in the other two formats: `"color_sub"`, `"color_df"`, and
#'    `"color_list"`.
#'
#' @examples
#' df <- data.frame(
#'    genotype=rep(c("WT", "GeneAKO", "GeneBKO"), c(4, 8, 8)),
#'    treatment=rep(rep(c("control", "treated"), each=2), 5),
#'    class=rep(c("WT", "KO"), c(4, 16)),
#'    time=c(rep("early", 4),
#'       rep(rep(c("early", "late"), each=4), 2)))
#' df$sample_group <- jamba::pasteByRow(df[,c("genotype", "treatment", "time")])
#' df$sample_name <- jamba::makeNames(df$sample_group);
#' df
#'
#' dfc <- design2colors(df,
#'    group_colnames="genotype",
#'    lightness_colnames="treatment",
#'    class_colnames="class")
#'
#' # same as above except assign colors to columns and some values
#' dfc <- design2colors(df,
#'    group_colnames="genotype",
#'    lightness_colnames="treatment",
#'    class_colnames="class",
#'    preset="dichromat",
#'    color_sub=c(KO="firebrick3",
#'       treatment="navy",
#'       class="cyan",
#'       time="dodgerblue"))
#'
#' # same as above except assign specific group colors
#' dfc <- design2colors(df,
#'    group_colnames="genotype",
#'    lightness_colnames="treatment",
#'    class_colnames="class",
#'    preset="dichromat",
#'    color_sub=c(
#'       WT_WT="gold",
#'       KO_GeneAKO="firebrick3",
#'       KO_GeneBKO="dodgerblue",
#'       treatment="navy",
#'       time="darkorchid4",
#'       class="burlywood"))
#'
#' dfc2 <- design2colors(df,
#'    group_colnames="genotype",
#'    lightness_colnames=c("time", "treatment"),
#'    class_colnames="class",
#'    preset="dichromat")
#'
#' dfc3 <- design2colors(df,
#'    group_colnames=c("genotype"),
#'    lightness_colnames=c("time", "treatment"),
#'    class_colnames="genotype",
#'    rotate_phase=TRUE,
#'    preset="dichromat")
#'
#' df1 <- df;
#' df2 <- subset(df, time %in% "early");
#' df12 <- rbind(df1, df2);
#' dfc12 <- design2colors(df12,
#'    group_colnames="genotype",
#'    lightness_colnames=c("time", "treatment"),
#'    class_colnames="-class",
#'    preset="dichromat")
#'
#' @export
design2colors <- function
(x,
   group_colnames=NULL,
   lightness_colnames=NULL,
   class_colnames=NULL,
   preset="dichromat",
   phase=1,
   rotate_phase=-1,
   class_pad=1,
   end_hue_pad=2,
   desat=c(0, 0.4),
   dex=c(2, 5),
   color_sub=NULL,
   plot_type=c("table",
      "list",
      "none"),
   return_type=c("list",
      "df",
      "vector"),
   verbose=FALSE,
   ...)
{
   #
   if (length(x) == 0) {
      return(NULL)
   }

   # Convert SummarizedExperiment to data.frame using colData(x)
   if ("SummarizedExperiment" %in% class(x)) {
      x <- data.frame(check.names=FALSE,
         SummarizedExperiment::colData(x));
   }

   # validate arguments
   plot_type <- match.arg(plot_type);
   return_type <- match.arg(return_type);
   if (length(class_pad) != 1 || class_pad < 0 ) {
      class_pad <- 1;
   }
   desat <- rep(desat,
      length.out=2);
   dex <- rep(dex,
      length.out=2);

   # handle each argument of colnames
   if (is.numeric(group_colnames)) {
      group_colnames <- jamba::rmNA(colnames(x)[group_colnames]);
   }
   if (length(lightness_colnames) > 0) {
      if (is.numeric(lightness_colnames)) {
         lightness_colnames <- jamba::rmNA(colnames(x)[lightness_colnames]);
      }
   }
   if (length(class_colnames) > 0) {
      if (is.numeric(class_colnames)) {
         class_colnames <- jamba::rmNA(colnames(x)[class_colnames]);
      }
   }

   # sort by class, group, lightness to make downstream steps consistent
   x_input <- x;
   # iterate each column and convert to factor if needed
   all_colnames1 <- gsub("^[-]", "",
      c(class_colnames,
         group_colnames,
         lightness_colnames));
   for (xcol in all_colnames1) {
      if (!is.factor(x[[xcol]])) {
         x[[xcol]] <- factor(x[[xcol]],
            levels=unique(x[[xcol]]));
      }
   }
   x <- jamba::mixedSortDF(x,
      byCols=c(class_colnames,
         group_colnames,
         lightness_colnames));
   class_colnames <- gsub("^[-]", "", class_colnames);
   group_colnames <- gsub("^[-]", "", group_colnames);
   lightness_colnames <- gsub("^[-]", "", lightness_colnames);
   all_colnames <- c(class_colnames,
      group_colnames,
      lightness_colnames);

   # convert character to factor to apply order each appears
   for (icol in all_colnames) {
      if (!is.factor(x[[icol]])) {
         x[[icol]] <- factor(x[[icol]],
            levels=unique(x[[icol]]));
      }
   }

   # generate output per class, group, lightness values
   xlist <- list(
      class=jamba::pasteByRowOrdered(x[,class_colnames, drop=FALSE]),
      group=jamba::pasteByRowOrdered(x[,group_colnames, drop=FALSE]),
      lightness=jamba::pasteByRowOrdered(x[,lightness_colnames, drop=FALSE])
   )
   xlist <- jamba::rmNULL(xlist, nullValue="");
   xdf <- data.frame(check.names=FALSE,
      xlist);

   # add class_group and class_group_lightness column values
   xdf$class_group <- jamba::pasteByRowOrdered(xdf[,c("class", "group"), drop=FALSE]);
   xdf$class_group_lightness <- jamba::pasteByRowOrdered(xdf[,c("class", "group", "lightness"), drop=FALSE]);

   # sort table by class, group, lightness
   xdf <- jamba::mixedSortDF(xdf,
      byCols=c("class",
         "group",
         "lightness",
         "class_group",
         "class_group_lightness"));
   if (verbose) {
      jamba::printDebug("design2colors(): ",
         "full class, group, lightness df:");
      print(xdf);
   }
   udf <- unique(xdf);
   if (verbose) {
      jamba::printDebug("design2colors(): ",
         "unique class, group, lightness df:");
      print(udf);
   }

   gdf <- unique(udf[,c("class", "group", "class_group"), drop=FALSE]);
   gdf$real <- 1;
   gdf0 <- head(gdf, 1);
   gdf0[1,] <- NA
   if (verbose) {
      jamba::printDebug("design2colors(): ",
         "unique class, group df:");
      print(gdf);
   }

   if (length(unique(gdf$class)) > 1) {
      ibreaks <- jamba::breaksByVector(gdf$class)$breakPoints;
      if (length(ibreaks) > 1) {
         for (k in rev(ibreaks)) {
            gdf_list <- c(
               list(head(gdf, k)),
               rep(list(gdf0), class_pad),
               list(tail(gdf, -k)));
            gdf <- jamba::rbindList(gdf_list);
         }
      }
   }
   if (verbose) {
      jamba::printDebug("design2colors(): ",
         "expanded unique class, group df:");
      print(gdf);
   }

   # assign color hue
   n <- nrow(gdf);
   offset <- 0;
   hue_seq <- head(seq(from=0 + offset,
      to=360 + offset,
      length.out=n + end_hue_pad), n);
   gdf$hue <- hue_seq;
   # adjust color hue
   gdf$hue2 <- colorjam::hw2h(hue_seq,
      preset="dichromat");
   gdf <- subset(gdf, !is.na(real))

   # assign colors
   class_group_hue1 <- jamba::nameVector(gdf[,c("hue", "class_group")])
   class_group_hue <- colorjam::hw2h(class_group_hue1,
      preset=preset);
   if (all(as.character(gdf$class_group) %in% names(color_sub))) {
      class_group_color <- color_sub[as.character(gdf$class_group)];
   } else {
      class_group_color <- colorjam::rainbowJam(n=length(class_group_hue),
         hues=class_group_hue,
         phase=phase,
         preset=preset,
         ...);
      names(class_group_color) <- gdf$class_group;
   }
   gdf$class_group_color <- class_group_color;
   if (verbose) {
      jamba::printDebug("design2colors(): ",
         "expanded unique class, group df, with colors:");
      print(gdf);
   }

   # optionally rotate phase
   phase <- phase + rotate_phase;

   udf$class_group_hue <- class_group_hue[as.character(udf$class_group)]
   udf$class_group_color <- class_group_color[as.character(udf$class_group)]
   xdf$class_group_color <- class_group_color[as.character(xdf$class_group)];
   class_group_lightness_color <- NULL;
   if (any(duplicated(udf$class_group_color))) {
      class_group_lightness_color <- jamba::color2gradient(udf$class_group_color,
         dex=dex[1]);
      names(class_group_lightness_color) <- udf$class_group_lightness;
      udf$class_group_lightness_color <- class_group_lightness_color;
      xdf$class_group_lightness_color <- class_group_lightness_color[as.character(xdf$class_group_lightness)];
   } else {
      class_group_lightness_color <- jamba::nameVector(
         udf$class_group_color,
         as.character(udf$class_group_lightness));
      udf$class_group_lightness_color <- class_group_lightness_color;
      xdf$class_group_lightness_color <- class_group_lightness_color[as.character(xdf$class_group_lightness)];
   }

   kcolnames <- intersect(
      c("class_group_lightness_color",
         "class_group_color"),
      colnames(udf));

   ###############################################################
   # assign group colors when cardinality is appropriate
   if (length(rownames(x)) == 0 ||
         any(!grepl("^[0-9]+$", rownames(x)))) {
      # if rownames are entirely integers, do not assign categorical colors
      iter_colnames <- jamba::nameVector(colnames(x));
   } else {
      iter_colnames <- jamba::nameVector(c("rownames", colnames(x)));
   }
   new_colors <- lapply(iter_colnames, function(icol) {
      if ("rownames" %in% icol) {
         ivalues <- data.frame(`rownames`=rownames(x));
      } else {
         ivalues <- x[,icol, drop=FALSE]
      }
      # jamba::printDebug("icol: ", icol);
      # jamba::printDebug("      ivalues: ", ivalues[,1]);
      for (kcolname in kcolnames) {
         # jamba::printDebug("      kcolname: ", kcolname);
         idf <- unique(data.frame(check.names=FALSE,
            ivalues,
            xdf[,kcolname, drop=FALSE]));
         if (any(duplicated(idf[[icol]]))) {
            next;
         } else {
            kcolors <- jamba::nameVector(
               as.character(idf[[2]]),
               as.character(idf[[1]]));
            # printDebug(kcolors);
            return(kcolors);
         }
      }
      return(NULL);
   });
   # create color vector
   new_colors_v <- unlist(unname(new_colors));

   ###############################################################
   # now generate categorical colors for remaining columns
   if (any(lengths(new_colors) == 0)) {
      add_colors_v1 <- NULL;
      add_colnames <- names(new_colors)[lengths(new_colors) == 0];

      ######################################
      # if color_sub matches colname,
      # use that color with gradient effect
      colname_colnames <- NULL;
      if (any(add_colnames %in% names(color_sub))) {
         colname_colnames <- intersect(add_colnames, names(color_sub));
         add_colnames <- setdiff(add_colnames, colname_colnames);
         colname_colors <- lapply(colname_colnames, function(icol){
            if (verbose > 1) {
               jamba::printDebug("design2colors(): ",
                  "colname_icol: ", icol);
            }
            if (is.factor(x_input[[icol]])) {
               if (verbose > 1) {
                  jamba::printDebug("design2colors(): ", c("   is.factor=", "TRUE"), sep="");
               }
               ivalues <- levels(x_input[[icol]]);
            } else {
               if (verbose > 1) {
                  jamba::printDebug("design2colors(): ", c("   is.factor=", "FALSE"), sep="");
               }
               ivalues <- unique(as.character(x_input[[icol]]));
            }
            icolors <- jamba::nameVector(
               jamba::color2gradient(color_sub[[icol]],
                  dex=dex[2],
                  n=length(ivalues)),
               ivalues);
            if (verbose > 1) {
               jamba::printDebug("design2colors(): ",
                  "ivalues: ", ivalues);
               jamba::printDebug("design2colors(): ",
                  "icolors:")
               jamba::printDebugI(icolors);
            }
            icolors;
         });
         add_colors_v1 <- unlist(unname(colname_colors));
         add_colors_v1 <- add_colors_v1[!duplicated(names(add_colors_v1))];
      }

      ########################################
      # all other colors are assigned  colors
      add_values <- unique(unlist(lapply(add_colnames, function(icol) {
         if ("rownames" %in% icol) {
            unique(rownames(x))
         } else {
            if (is.factor(x_input[[icol]])) {
               levels(x_input[[icol]])
            } else {
               unique(as.character(x_input[[icol]]))
            }
         }
      })));
      # reuse color assignments if present
      if (any(c(add_values, names(add_colors_v1)) %in% names(color_sub))) {
         add_values_1 <- intersect(
            c(add_values,
               names(add_colors_v1)),
            names(color_sub));
         add_values <- setdiff(add_values, add_values_1);
         add_colors_v1 <- c(
            color_sub[add_values_1],
            add_colors_v1)
         add_colors_v1 <- add_colors_v1[!duplicated(names(add_colors_v1))];
      }

      # sometimes a factor is already assigned a color
      if (any(add_values %in% names(new_colors_v))) {
         add_values <- setdiff(add_values,
            names(new_colors_v))
      }
      add_n <- length(add_values);

      if (add_n > 0) {
         # determine "improved" starting hue using mean of first two hues
         if (nrow(gdf) > 1) {
            start_hue1 <- mean(gdf[1:2, "hue"])
         } else {
            start_hue1 <- 0;
         }
         offset <- 0;
         add_hue_seq <- head(seq(from=0 + start_hue1,
            to=360 + start_hue1,
            length.out=add_n + end_hue_pad), add_n);
         add_hue <- colorjam::hw2h(add_hue_seq,
            preset=preset);
         add_hue1 <- head(as.vector(
            matrix(ncol=2,
               byrow=TRUE,
               add_hue)),
            add_n);

         # assign categorical colors to unique values
         add_m <- head(as.vector(
            matrix(ncol=2,
               byrow=TRUE,
               add_values)),
            add_n);
         add_m <- add_values;
         add_colors_v <- jamba::nameVector(
            colorjam::rainbowJam(add_n,
               phase=phase,
               hues=add_hue,
               preset=preset,
               ...),
            add_m);
      } else {
         add_colors_v <- NULL;
      }
      add_colors_v <- c(add_colors_v,
         add_colors_v1);
      add_colors_v <- add_colors_v[!duplicated(names(add_colors_v))];
      add_colors_v <- jamba::nameVector(
         colorspace::desaturate(add_colors_v,
            amount=desat[2]),
         names(add_colors_v));

      # optionally rotate phase
      phase <- phase + rotate_phase;

      all_add_colors_v <- c(add_colors_v,
         new_colors_v);

      add_colors <- lapply(jamba::nameVector(c(add_colnames, colname_colnames)), function(icol) {
         all_add_colors_v[unique(as.character(x[[icol]]))]
      })
   } else {
      add_colors <- NULL;
   }
   # end filling in new_colors
   ###############################

   all_colors_list1 <- jamba::rmNULL(
      c(
         c(jamba::rmNULL(new_colors), add_colors)[colnames(x)],
         list(class_group_color=class_group_color,
            class_group_lightness_color=class_group_lightness_color)));
   all_colors_v <- unlist(unname(all_colors_list1));
   if (any(duplicated(names(all_colors_v)))) {
      all_colors_v <- all_colors_v[!duplicated(names(all_colors_v))];
   }
   all_colors_list <- lapply(all_colors_list1, function(i){
      all_colors_v[names(i)];
   })
   if (verbose > 1) {
      jamba::printDebug("design2colors(): ",
         "all_colors_list1: ");
      print(all_colors_list1);
      jamba::printDebug("design2colors(): ",
         "all_colors_v: ");
      jamba::printDebugI(all_colors_v);
      jamba::printDebug("design2colors(): ",
         "all_colors_list: ");
      print(all_colors_list);
   }

   if (!desat[1] == 0) {
      all_colors_list <- lapply(all_colors_list, function(i){
         jamba::nameVector(
            colorspace::desaturate(col=i,
               amount=desat[1]),
            names(i))
      })
   }

   # one option is to display the full list of colors
   if ("list" %in% plot_type) {
      jamba::showColors(all_colors_list)
   }

   # another option is to display the input data.frame colorized
   x_colors_list <- lapply(jamba::nameVector(colnames(x_input)), function(i){
         all_colors_list[[i]][as.character(x_input[[i]])]
   });
   x_colors <- as.data.frame(
      lapply(jamba::nameVector(colnames(x_input)), function(i){
         all_colors_list[[i]][as.character(x_input[[i]])]
      }));
   if ("table" %in% plot_type) {
      opar <- par(no.readonly=TRUE);
      jamba::adjustAxisLabelMargins(
         rownames(x_colors), 2)
      jamba::imageByColors(x_colors,
         cellnote=x_input,
         flip="y",
         cexCellnote=0.7)
      par(opar);
   }

   if ("df" %in% return_type) {
      attr(x_colors, "color_list") <- all_colors_list;
      attr(x_colors, "color_sub") <- all_colors_v;
      attr(x_colors, "df") <- x_input;
      return(x_colors)
   }
   if ("vector" %in% return_type) {
      attr(all_colors_v, "color_list") <- all_colors_list;
      attr(all_colors_v, "df") <- x_input;
      attr(all_colors_v, "color_df") <- x_colors;
      return(all_colors_v);
      #return(unlist(unname(all_colors_list)));
   }

   attr(all_colors_list, "df") <- x_input;
   attr(all_colors_list, "color_df") <- x_colors;
   attr(all_colors_list, "color_sub") <- all_colors_v;
   return(all_colors_list)
}


#' Determine cardinality between two vectors
#'
#' Determine cardinality between two vectors
#'
#' @examples
#' a <- letters[c(1, 1, 2, 2, 2, 3)];
#' b <- LETTERS[c(1, 2, 3, 4, 5, 6)];
#' d <- LETTERS[c(1, 2, 2, 1, 1, 1)];
#'
#' cardinality(a, b)
#' cardinality(b, a)
#'
#' ab <- data.frame(a, b)
#' ad <- data.frame(a, d)
#' cardinality(ab, ad)
#'
#' abt <- tibble::tibble(a, b);
#' adt <- tibble::tibble(a, d);
#' cardinality(abt, adt)
#'
#' abm <- as.matrix(abt);
#' adm <- as.matrix(adt);
#' cardinality(abm, adm)
#'
#' cardinality(d, adm, verbose=TRUE)
#'
#' cardinality(d, adm, verbose=2)
#'
#' @export
cardinality <- function
(x,
   y=NULL,
   verbose=FALSE,
   ...)
{
   #
   df_classes <- c("matrix",
      "data.frame",
      "DataFrame",
      "DFrame",
      "tbl");
   if (length(y) > 0 && any(df_classes %in% class(y))) {
      if (verbose) {
         jamba::printDebug("cardinality(): ",
            c("Converting y to vector with ", "pasteByRow()"),
            sep="");
      }
      y <- jamba::pasteByRow(y);
   }
   if (length(x) > 0) {
      if (any(df_classes %in% class(x))) {
         if (length(y) > 0) {
            if (verbose) {
               jamba::printDebug("cardinality(): ",
                  c("Converting x to vector with ", "pasteByRow()"),
                  sep="");
            }
            x <- jamba::pasteByRow(x);
            x <- data.frame(x=x,
               y=y);
         } else {
            if (!ncol(x) == 2) {
               stop("When x is supplied alone, it must have two columns.");
            }
            x <- data.frame(x=x[,1],
               y=x[,2]);
         }
      } else {
         if (length(y) == 0) {
            stop("When x is a vector, y must be supplied.");
         }
         x <- data.frame(x=x,
            y=y);
      }
   }

   # x is a data.frame with two columns
   x_uniq <- unique(x);

   if (verbose > 1) {
      jamba::printDebug("cardinality(): ",
         "head(x, 20):")
      print(head(x, 20));
      jamba::printDebug("cardinality(): ",
         "head(unique(x), 20):")
      print(head(x_uniq, 20));
   }

   x_tc <- jamba::tcount(x_uniq[,1]);
   y_tc <- jamba::tcount(x_uniq[,2]);
   c(`from`=max(x_tc),
      `to`=max(y_tc))
}

#' Convert data.frame to numeric color substitution
#'
#' Convert data.frame to numeric color substitution
#'
#' @export
df_to_numcolors <- function
(df,
   colramp=c("Reds"),
   colramp_divergent=c("RdBu_r"),
   ...)
{
   # find numeric columns
   numcols <- which(unname(jamba::sclass(df)) %in% c("integer", "numeric", "float"))
   colramp <- rep(colramp, length.out=length(numcols));
   colramp_divergent <- rep(colramp_divergent, length.out=length(numcols));

   # calculate max value per column
   nummaxs <- sapply(numcols, function(i){
      max(abs(df[[i]]), na.rm=TRUE)
   });

   # iterate each column and populate color vector
   colsub <- list();
   for (i in rev(seq_along(numcols))) {
      vals <- df[[numcols[i]]]
      if (nummaxs[i] > 10) {
         vals <- round(vals);
         df[[numcols[i]]] <- vals;
      } else {
         vals <- round(vals * 100) / 100
         df[[numcols[i]]] <- vals;
      }
      vals <- unique(rmNA(vals));
      if (min(vals) < 0) {
         col1 <- colorjam::col_div_xf(max(abs(vals)),
            colramp=colramp_divergent)
         col1sub <- jamba::nameVector(col1(vals), vals)
      } else {
         col1 <- circlize::colorRamp2(
            colors=jamba::getColorRamp("Reds", trimRamp=c(2, 2), n=5),
            breaks=seq(from=min(vals)-0.1, to=max(vals)+1, length.out=5))
         col1sub <- jamba::nameVector(col1(vals), vals)
      }
      colsub[as.character(vals)] <- col1sub;
   }
   return(list(
      df=df,
      color_sub=colsub));
}
