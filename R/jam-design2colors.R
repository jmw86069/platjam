
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
#' @param force_consistent_colors `logical` indicating whether to force
#'    color substitutions across multiple columns, when those columns
#'    share one or more of the same values. Note: This scenario is most
#'    likely to occur when using `color_sub` to assign colors to a
#'    specific column in `colnames(x)`, and where that column may contain
#'    one or more values already assigned a color earlier in the process.
#'    For example: class/group/lightness defines colors for these columns,
#'    then all columns are checked for cardinality with these color
#'    assignments. Columns not appearing in `color_sub` will be colorized
#'    when their cardinality is compatible with group colors, otherwise
#'    specific values may be assigned colors via `color_sub`, then all
#'    remaining values are assigned categorical colors. This process defines
#'    colors for each column. The last step reassigns colors consistently
#'    using the first value-color association that appears in the list,
#'    to make sure all assignments are consistent. This last step is subject
#'    of the argument `force_consistent_colors`. In either case, the
#'    output `color_sub` will only have one color assigned per value.
#'    * Default `TRUE`: When colors are being assigned to column values
#'    `color_sub`, if the value had been assigned a color in a previous
#'    column, for example by `group_colnames` color assignment, then the
#'    first assignment is applied to all subsequent assignments.
#'    * Optional `FALSE`: Color assignments are re-used where applicable,
#'    except when overridden by `color_sub` for a particular column. In
#'    that case, color assignments are maintained for each specific column.
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
#' df$age <- sample(40:80, size=nrow(df));
#' df
#'
#' dfc <- design2colors(df,
#'    group_colnames="genotype",
#'    lightness_colnames="treatment",
#'    class_colnames="class",
#'    color_sub=c(age="dodgerblue"))
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
#'       WT="gold",
#'       KO="purple3",
#'       GeneAKO="firebrick3",
#'       GeneBKO="dodgerblue",
#'       treatment="navy",
#'       time="darkorchid4"))
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
#'    class_colnames="class",
#'    preset="dichromat",
#'    color_sub=c(
#'       treatment="steelblue",
#'       time="dodgerblue"
#'    ))
#'
#'
#' @export
design2colors <- function
(x,
   group_colnames=NULL,
   lightness_colnames=NULL,
   class_colnames=NULL,
   ignore_colnames=NULL,
   preset="dichromat",
   phase=1,
   rotate_phase=-1,
   class_pad=1,
   end_hue_pad=2,
   hue_offset=0,
   desat=c(0, 0.4),
   dex=c(2, 5),
   color_sub=NULL,
   force_consistent_colors=TRUE,
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

   # optionally ignore some colnames
   if (any(ignore_colnames %in% colnames(x))) {
      x <- x[,setdiff(colnames(x), ignore_colnames), drop=FALSE];
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

   # check cardinality of class and group
   if (length(class_colnames) > 0 && length(group_colnames)) {
      class_card <- cardinality(x[,class_colnames, drop=FALSE],
         x[,group_colnames, drop=FALSE]);
      if (class_card["from"] != 1) {
         stop("class_colnames must have 1-to-X cardinality with group_colnames.");
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
   # xdf$class_group <- jamba::pasteByRowOrdered(xdf[,c("class", "group"), drop=FALSE]);
   # xdf$class_group_lightness <- jamba::pasteByRowOrdered(xdf[,c("class", "group", "lightness"), drop=FALSE]);
   xdf$class_group <- jamba::pasteByRowOrdered(xdf[,c("group"), drop=FALSE]);
   xdf$class_group_lightness <- jamba::pasteByRowOrdered(xdf[,c("group", "lightness"), drop=FALSE]);

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
   rownames(gdf0) <- "class_pad";
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
   #hue_offset <- 0;
   hue_seq <- head(seq(from=0 + hue_offset,
      to=360 + hue_offset,
      length.out=n + end_hue_pad), n);
   gdf$hue <- hue_seq;
   # adjust color hue
   #gdf$hue2 <- colorjam::hw2h(hue_seq,
   #   preset="dichromat");
   gdf <- subset(gdf, !is.na(real))

   # assign colors
   class_group_hue1 <- jamba::nameVector(gdf[,c("hue", "class_group")])
   class_group_hue <- colorjam::hw2h(class_group_hue1,
      preset=preset);
   if (all(as.character(gdf$class_group) %in% names(color_sub))) {
      class_group_color <- color_sub[as.character(gdf$class_group)];
      gdf$hue <- colorjam::h2hw(jamba::col2hcl(class_group_color)["H",],
         preset=preset);
   } else {
      class_group_color <- colorjam::rainbowJam(n=length(class_group_hue),
         hues=class_group_hue,
         phase=phase,
         preset=preset,
         ...);
      names(class_group_color) <- gdf$class_group;
   }
   gdf$class_group_color <- class_group_color;

   # get mean color hue per class
   class_hues <- NULL;
   if (length(class_colnames) > 0) {
      # calculate mean hue per class
      class_hues <- sapply(split(gdf$hue, gdf$class), function(ihue){
         if (length(ihue) %% 2 > 0) {
            ihue[ceiling(length(ihue)/2)] %% 360;
         } else {
            colorjam::h2hw(
               mean_hue(colorjam::hw2h(ihue,
                  preset=preset)),
               preset=preset)
         }
      })
      class_hue <- colorjam::hw2h(class_hues,
         preset=preset);
      # class_hue <- class_hues;
      if (all(as.character(gdf$class) %in% names(color_sub))) {
         class_color <- color_sub[as.character(gdf$class)];
      } else {
         class_color <- colorjam::rainbowJam(n=length(class_hue),
            hues=class_hue,
            phase=phase + 3,
            preset=preset,
            ...);
         names(class_color) <- names(class_hues);
         color_sub[names(class_color)] <- class_color;
      }
   }

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

   ################################################################
   # assign group colors when cardinality is appropriate
   # Note: use all columns including class, group, lightness
   # which allows various combinations to be detected and assigned.
   if (length(rownames(x)) == 0 ||
         any(!grepl("^[0-9]+$", rownames(x)))) {
      # if rownames are entirely integers, do not assign categorical colors
      iter_colnames <- jamba::nameVector(colnames(x));
   } else {
      iter_colnames <- jamba::nameVector(c("rownames", colnames(x)));
   }
   new_colors <- lapply(iter_colnames, function(icol) {
      # new in version 0.0.52.900: do not assign colors when color_sub is defined
      if (length(color_sub) > 0 && icol %in% names(color_sub)) {
         return(NULL);
      }
      # skip numeric columns which will be assigned gradient colors
      if (is.numeric(x[[icol]])) {
         return(NULL);
      }
      if ("rownames" %in% icol) {
         ivalues <- data.frame(`rownames`=rownames(x));
      } else {
         ivalues <- x[,icol, drop=FALSE]
      }
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
   new_color_list <- c(new_colors);
   color_sub_1 <- color_sub;
   color_sub <- c(color_sub_1,
      unlist(unname(new_color_list)))

   ############################################
   # now generate colors for remaining columns
   add_color_functions <- NULL;
   colname_colors <- NULL;
   # any empty element in new_colors lacks color assignment
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

         # iterate remaining colnames and assign colors
         colname_colors <- lapply(jamba::nameVector(colname_colnames), function(icol){
            if (verbose > 1) {
               jamba::printDebug("design2colors(): ",
                  "colname_icol: ", icol);
            }
            if (is.factor(x_input[[icol]])) {
               if (verbose > 1) {
                  jamba::printDebug("design2colors(): ", c("   is.factor=", "TRUE"), sep="");
               }
               ivalues <- levels(x_input[[icol]]);
            } else if (is.numeric(x_input[[icol]])) {
               # numeric columns will receive gradient color function
               if (verbose > 1) {
                  jamba::printDebug("design2colors(): ", c("   is.numeric=", "TRUE"), sep="");
               }
               icolors <- assign_numeric_colors(x=x_input[[icol]],
                  restrict_pretty_range=FALSE,
                  color=color_sub[[icol]]);
               return(icolors);
            } else {
               if (verbose > 1) {
                  jamba::printDebug("design2colors(): ", c("   is.factor=", "FALSE"), sep="");
               }
               ivalues <- unique(as.character(x_input[[icol]]));
            }
            # use ivalues to define colors
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
         if (verbose > 1) {
            jamba::printDebug("design2colors(): ",
               "colname_colors: ");
            print_color_list(colname_colors);
         }
         new_color_list[names(colname_colors)] <- colname_colors;
         color_sub <- c(color_sub_1,
            unlist(unname(new_color_list[!jamba::sclass(new_color_list) %in% "function"])))

         # extract new color assignments so they can be consistent
         # with values in other columns if needed.
         # Note: This step ignores color functions for numeric columns.
         # is_color_fn <- sapply(colname_colors, is.function);
         # add_colors_v1 <- unlist(unname(colname_colors[!is_color_fn]));
         # add_colors_v1 <- add_colors_v1[!duplicated(names(add_colors_v1))];
         # if (any(is_color_fn)) {
         #    add_color_functions <- colname_colors[is_color_fn];
         #    names(add_color_functions) <- colname_colnames[is_color_fn];
         #    if (length(add_color_functions) > 0 && verbose > 1) {
         #       jamba::printDebug("design2colors(): ",
         #          "add_color_functions: ", names(add_color_functions));
         #    }
         # }
      }

      ###################################################
      # all other values are assigned categorical colors
      if (verbose > 1) {
         jamba::printDebug("design2colors(): ",
            "add_colnames: ", add_colnames);
      }
      # assemble all remaining column values for color assignment
      add_numeric_colnames <- add_colnames[sapply(add_colnames, function(icol){
         is.numeric(x_input[[icol]])})];
      # assemble all unique values that need colors,
      # using column name itself for numeric columns
      add_values <- unique(unlist(
         lapply(add_colnames, function(icol) {
            if ("rownames" %in% icol) {
               unique(rownames(x))
            } else {
               if (is.factor(x_input[[icol]])) {
                  levels(x_input[[icol]])
               } else if (is.numeric(x_input[[icol]])) {
                  # for numeric columns assign a color to the colname
                  icol;
               } else {
                  unique(as.character(x_input[[icol]]))
               }
            }
         })
      ));
      # reuse color assignments if present
      # if (any(c(add_values, names(add_colors_v1)) %in% names(color_sub))) {
      # Note: do not reassign names(add_colors_v1) via names(color_sub).
      # Give preference to new_colors_v, defined by class/group/lightness
      # remove add_values that already have assignments in color_sub
      if (any(add_values %in% names(color_sub))) {
         add_values_1 <- intersect(
            add_values,
            names(color_sub));
         add_values <- setdiff(add_values, add_values_1);
         # add_match <- match(add_values_1,
         #    c(names(new_colors_v),
         #       names(color_sub)));
         # add_colors_v1 <- c(new_colors_v,
         #    color_sub)[add_match];
      }
      add_n <- length(add_values);

      # check if any values remain to be assigned new colors
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
         # optionally "rotate" the assignment by column to avoid
         # adjacent colors being too similar
         # add_m <- as.vector(
         #    matrix(ncol=2,
         #       byrow=TRUE,
         #       add_values));
         # add_m <- add_m[!duplicated(add_m)];
         add_m <- add_values;
         add_colors_v <- jamba::nameVector(
            colorjam::rainbowJam(add_n,
               phase=phase,
               hues=add_hue,
               preset=preset,
               ...),
            add_m);
         # extend color_sub with newly assigned colors
         color_sub <- c(color_sub,
            add_colors_v);

         # now re-apply these color_sub to each add_colnames
         add_colname_colors <- lapply(jamba::nameVector(add_colnames), function(icol){
            if ("rownames" %in% icol) {
               icolors <- color_sub[unique(rownames(x))]
            } else {
               if (is.factor(x_input[[icol]])) {
                  icolors <- color_sub[levels(x_input[[icol]])]
               } else if (is.numeric(x_input[[icol]])) {
                  icolors <- assign_numeric_colors(x=x_input[[icol]],
                     restrict_pretty_range=FALSE,
                     color=color_sub[[icol]]);
               } else {
                  icolors <- color_sub[unique(as.character(x_input[[icol]]))]
               }
            }
            icolors;
         })
         if (verbose > 1) {
            jamba::printDebug("design2colors(): ",
               "add_colname_colors: ");
            print(sdim(add_colname_colors));
            print(add_colname_colors);
            print_color_list(add_colname_colors);
         }
         new_color_list[names(add_colname_colors)] <- add_colname_colors;
         color_sub <- c(color_sub,
            unlist(unname(
               add_colname_colors[!jamba::sclass(add_colname_colors) %in% "function"])))
      } else {
         # no new colors need to be assigned
         # add_colors_v <- NULL;
      }
      # optionally rotate phase
      phase <- phase + rotate_phase;

   } else {
      add_colors <- NULL;
      add_color_functions <- NULL;
   }
   # end filling in new_colors
   ###############################

   ###############################
   # Refresh the list of colors
   all_colors_list1 <- jamba::rmNULL(
      c(new_color_list[colnames(x_input)],
         list(
            class_group_color=class_group_color,
            class_group_lightness_color=class_group_lightness_color)));
   is_color_fn <- sapply(all_colors_list1, is.function);

   # obtain vector of value-color associations except for columns
   # that use color functions
   all_colors_v <- unlist(unname(all_colors_list1[!is_color_fn]));
   # remove any duplicated names, where the same value
   # appeared in multiple columns. This step retains only the
   # first color assignment per value.
   if (any(duplicated(names(all_colors_v)))) {
      all_colors_v <- all_colors_v[!duplicated(names(all_colors_v))];
   }

   # optionally force consistent color assignment by column
   if (force_consistent_colors) {
      all_colors_list <- lapply(all_colors_list1, function(i){
         if (is.function(i)) {
            i
         } else {
            all_colors_v[names(i)];
         }
      })
   } else {
      all_colors_list <- all_colors_list1;
   }

   if (verbose > 1) {
      # jamba::printDebug("design2colors(): ",
      #    "add_colors: ");
      # print(add_colors);
      jamba::printDebug("design2colors(): ",
         "all_colors_list1: ");
      print(sdim(all_colors_list1));
      print_color_list(all_colors_list1);
      jamba::printDebug("design2colors(): ",
         "all_colors_v: ");
      print_color_list(all_colors_v);
      jamba::printDebug("design2colors(): ",
         "all_colors_list: ");
      print_color_list(all_colors_list);
   }

   if (!desat[1] == 0) {
      all_colors_list <- lapply(all_colors_list, function(i){
         if (is.function(i)) {
            i;
         } else {
            jamba::nameVector(
               colorspace::desaturate(col=i,
                  amount=desat[1]),
               names(i))
         }
      })
   }

   # one option is to display the full list of colors
   if ("list" %in% plot_type) {
      all_colors_list_use <- print_color_list(all_colors_list,
         do_print=FALSE);
      jamba::showColors(all_colors_list_use)
   }

   # another option is to display the input data.frame colorized
   x_colors_list <- lapply(jamba::nameVector(colnames(x_input)), function(i){
      if (is.function(all_colors_list[[i]])) {
         jamba::nameVector(all_colors_list[[i]](x_input[[i]]),
            round(x_input[[i]],
               digits=3));
      } else {
         all_colors_list[[i]][as.character(x_input[[i]])]
      }
   });
   x_colors_list <- (
      lapply(jamba::nameVector(colnames(x_input)), function(i){
         if (is.function(all_colors_list[[i]])) {
            jamba::nameVector(all_colors_list[[i]](x_input[[i]]),
               round(x_input[[i]],
                  digits=3));
         } else {
            all_colors_list[[i]][as.character(x_input[[i]])]
         }
      }));
   x_colors <- as.data.frame(x_colors_list);
   if ("table" %in% plot_type) {
      opar <- par(no.readonly=TRUE);
      jamba::adjustAxisLabelMargins(
         x=rownames(x_colors),
         margin=2)
      jamba::adjustAxisLabelMargins(
         x=colnames(x_colors),
         margin=1)
      jamba::imageByColors(x_colors,
         cellnote=if(nrow(x_colors) < 1000){ x_input} else {NULL},
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

   x_tc <- jamba::tcount(x_uniq[,2]);
   y_tc <- jamba::tcount(x_uniq[,1]);
   c(`from`=max(x_tc),
      `to`=max(y_tc))
}

#' Convert data.frame to numeric color substitution
#'
#' Convert data.frame to numeric color substitution
#'
#' @param df `data.frame` where columns with numeric values will be
#'    colorized.
#' @param colramp,colramp_divergent `character` name of color ramp for
#'    linear and divergent color gradients, respectively.
#' @param trimRamp `integer` vector length=2, passed to
#'    `jamba::getColorRamp()`, to trim the edge colors from a linear
#'    color gradient, thus avoiding extreme colors.
#' @param ... additional arguments are ignored.
#'
#' @export
df_to_numcolors <- function
(df,
 colramp=c("Reds"),
 colramp_divergent=c("RdBu_r"),
 trimRamp=c(1, 2),
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
      if (nummaxs[i] > 100) {
         vals <- round(vals);
         df[[numcols[i]]] <- vals;
      } else if (nummaxs[i] > 10) {
         vals <- round(vals * 10) / 10;
         df[[numcols[i]]] <- vals;
      } else {
         vals <- round(vals * 100) / 100
         df[[numcols[i]]] <- vals;
      }
      vals <- unique(rmNA(vals));
      if (min(vals) < 0) {
         col1 <- colorjam::col_div_xf(max(abs(vals)),
            colramp=colramp_divergent[[i]])
         col1sub <- jamba::nameVector(col1(vals), vals)
      } else {
         col1 <- circlize::colorRamp2(
            colors=jamba::getColorRamp(colramp[[i]],
               trimRamp=trimRamp,
               n=5),
            breaks=seq(from=min(vals), to=max(vals) + 1, length.out=5))
         col1sub <- jamba::nameVector(col1(vals), vals)
      }
      colsub[as.character(vals)] <- col1sub;
   }
   return(list(
      df=df,
      color_sub=colsub));
}
