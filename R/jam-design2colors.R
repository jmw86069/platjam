
#' Convert experiment design into categorical colors
#'
#' Convert experiment design into categorical colors
#'
#' @param x `data.frame` with columns to be colorized
#' @param group_colnames `character` or `intger` vector indicating
#'    which `colnames(x)` to use, in order, for group color assignment.
#' @param lightness_colnames `character` or `intger` vector indicating
#'    which `colnames(x)` to use, in order, for group lightness gradient.
#' @param class_colnames `character` or `intger` vector indicating
#'    higher-level grouping of `group_colnames`
#' @param ... additional arguments are passed to downstream functions.
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
#'    class_colnames="class",
#'    preset="rgb")
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
   desat=0,
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

   # validate arguments
   plot_type <- match.arg(plot_type);
   return_type <- match.arg(return_type);
   if (length(class_pad) != 1 || class_pad < 0 ) {
      class_pad <- 1;
   }
   desat <- rep(desat,
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

   xlist <- list(
      class=jamba::pasteByRowOrdered(x[,class_colnames, drop=FALSE]),
      group=jamba::pasteByRowOrdered(x[,group_colnames, drop=FALSE]),
      lightness=jamba::pasteByRowOrdered(x[,lightness_colnames, drop=FALSE])
   )
   xlist <- jamba::rmNULL(xlist, nullValue="");
   xdf <- data.frame(check.names=FALSE,
      xlist);
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
   class_group_color <- colorjam::rainbowJam(n=length(class_group_hue),
      hues=class_group_hue,
      phase=phase,
      preset=preset,
      ...);
   names(class_group_color) <- gdf$class_group;
   gdf$class_group_color <- class_group_color;

   # optionally rotate phase
   phase <- phase + rotate_phase;

   udf$class_group_hue <- class_group_hue[as.character(udf$class_group)]
   udf$class_group_color <- class_group_color[as.character(udf$class_group)]
   xdf$class_group_color <- class_group_color[as.character(xdf$class_group)];
   class_group_lightness_color <- NULL;
   if (any(duplicated(udf$class_group_color))) {
      class_group_lightness_color <- jamba::color2gradient(udf$class_group_color);
      names(class_group_lightness_color) <- udf$class_group_lightness;
      udf$class_group_lightness_color <- class_group_lightness_color;
      xdf$class_group_lightness_color <- class_group_lightness_color[as.character(xdf$class_group_lightness)];
   }

   kcolnames <- intersect(
      c("class_group_lightness_color",
         "class_group_color"),
      colnames(udf));

   # assign group colors when cardinality is appropriate
   if (length(rownames(x)) == 0 ||
         !any(is.na(as.numeric(rownames(x))))) {
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
      # jamba::printDebug("      ivalues: ", ivalues);
      for (kcolname in kcolnames) {
         # jamba::printDebug("      kcolname: ", kcolname);
         idf <- unique(data.frame(check.names=FALSE,
            ivalues,
            xdf[,kcolname, drop=FALSE]));
         if (any(duplicated(idf[[icol]]))) {
            next;
         } else {
            kcolors <- jamba::nameVector(idf[[2]],
               as.character(idf[[1]]));
            return(kcolors);
         }
      }
      return(NULL);
   });
   # create color vector
   new_colors_v <- unlist(unname(new_colors));
   # jamba::printDebug("new_colors_v:");
   # jamba::printDebug(new_colors_v);

   # now generate categorical colors for remaining columns
   if (any(lengths(new_colors) == 0)) {
      add_colnames <- names(new_colors)[lengths(new_colors) == 0];
      add_values <- unique(unlist(lapply(add_colnames, function(icol) {
         if ("rownames" %in% icol) {
            unique(rownames(x))
         } else {
            unique(x[[icol]])
         }
      })));
      # sometimes a factor is already assigned a color
      if (any(add_values %in% names(new_colors_v))) {
         add_values <- setdiff(add_values,
            names(new_colors_v))
      }
      add_n <- length(add_values);

      # determine "improved" starting hue using mean of first two hues
      if (nrow(gdf) > 1) {
         start_hue1 <- mean(gdf[1:2, "hue"])
      } else {
         start_hue1 <- 0;
      }
      # jamba::printDebug("start_hue1: ", start_hue1);
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
      add_colors_v <- jamba::nameVector(
         colorspace::desaturate(add_colors_v,
            amount=desat[2]),
         add_m);
      # jamba::printDebug("add_colors_v:");
      # jamba::printDebugI(add_colors_v);

      # optionally rotate phase
      phase <- phase + rotate_phase;

      all_add_colors_v <- c(add_colors_v,
         new_colors_v);

      add_colors <- lapply(jamba::nameVector(add_colnames), function(icol) {
         all_add_colors_v[as.character(unique(x[[icol]]))]
      })
   }
   all_colors_list <- jamba::rmNULL(c(
      c(jamba::rmNULL(new_colors), add_colors)[colnames(x)],
      list(class_group_color=class_group_color,
         class_group_lightness_color=class_group_lightness_color)));

   if (!desat == 0) {
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
   x_colors <- as.data.frame(
      lapply(jamba::nameVector(colnames(x)), function(i){
         all_colors_list[[i]][x[[i]]]
      }));
   if ("table" %in% plot_type) {
      opar <- par(no.readonly=TRUE);
      jamba::adjustAxisLabelMargins(
         rownames(x_colors), 2)
      jamba::imageByColors(x_colors,
         cellnote=x,
         flip="y",
         cexCellnote=0.7)
      par(opar);
   }

   if ("df" %in% return_type) {
      return(x_colors)
   }
   if ("vector" %in% return_type) {
      return(unlist(unname(all_colors_list)));
   }

   return(all_colors_list)
}
