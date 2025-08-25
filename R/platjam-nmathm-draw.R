
# function to check whether a character string is a color, color ramp or function

# function to draw the output from nmatlist2heatmaps()

#' Draw heatmaps using output from nmatlist2heatmaps()
#'
#' Draw heatmaps using output from nmatlist2heatmaps()
#'
#' @param nmathm `list` output from `nmatlist2heatmaps()`.
#' @param hm_nrow `integer` to define the number of rows to display heatmaps,
#'    currently experimental. In general, multi-row heatmap output is not
#'    easy to define as general-use visualization.
#'    Should color keys be shared by all rows?
#'    Should heatmaps break at exactly half the total number of heatmaps?
#'    Should row annotations be displayed on each row of heatmaps? And only
#'    display the color key for the first row of heatmaps?
#'    It is complicated, and not yet well-explored.
#' @param do_plot `logical` default TRUE, whether to draw the heatmap, or
#'    just prepare the object(s) to be drawn.
#'    Currently experimental, output formats are not stable yet.
#'    This option is mainly a placeholder to be able to run this function
#'    and produce *something* without actually rendering.
#' @param title,caption `character` with optional overall title and caption.
#'    By default it uses `nmathm$metadata$title` and `nmathm$metadata$caption`
#'    as already defined, but these options permit custom values.
#' @param ... additional arguments are passed to internal functions.
#'
#' @param
#'
#' @export
draw_nmathm <- function
(nmathm,
 hm_nrow=1,
 do_plot=TRUE,
 title=NULL,
 caption=NULL,
 ...)
{
   #

   HM_temp <- nmathm$draw$HM_temp;
   ht_gap <- nmathm$draw$ht_gap;
   main_heatmap <- nmathm$draw$main_heatmap;
   main_heatmap_temp <- main_heatmap;
   ht_gap <- nmathm$fn_params$ht_gap;
   title_gp <- nmathm$fn_params$title_gp;

   use_annotation_legend_list <- nmathm$fn_params$use_annotation_legend_list

   # grab the title
   if ("metadata" %in% names(nmathm)) {
      #
      if (length(title) == 0) {
         title <- nmathm$metadata$title;
      }
      if (length(caption) == 0) {
         caption <- nmathm$metadata$caption;
      }
   }
   if (do_plot &&
         (length(title) > 0 || length(caption) > 0)) {
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "Preparing ComplexHeatmap::draw(HeatmapList)");
      }
      HM_drawn <- ComplexHeatmap::draw(HM_temp,
         column_title=title,
         column_title_gp=title_gp,
         ht_gap=ht_gap,
         adjust_annotation_extension=TRUE,
         annotation_legend_list=use_annotation_legend_list,
         main_heatmap=main_heatmap_temp,
         merge_legends=TRUE,
         padding=padding)
      #
   } else if (do_plot) {
      HM_drawn <- ComplexHeatmap::draw(HM_temp,
         ht_gap=ht_gap,
         adjust_annotation_extension=TRUE,
         annotation_legend_list=use_annotation_legend_list,
         main_heatmap=main_heatmap_temp,
         merge_legends=TRUE,
         padding=padding)
   }
   ######################################
   ## Layout heatmap panels
   ##
   HM_drawn <- NULL;
   if (length(hm_nrow) > 0 &&
         hm_nrow > 1) {
      ######################################
      ## Optional multi-row layout
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "Applying multi-row layout.");
      }
      hm_split <- rep(
         rep(
            seq_len(hm_nrow),
            each=ceiling(length(nmatlist) / hm_nrow)),
         length.out=length(nmatlist));
      EH_l3 <- split(EH_l, hm_split);
      HM_temp <- NULL;
      main_heatmap_temp <- head(main_heatmap, 1) +
         (length(partition) > 0) +
         (length(AHM) > 0);
      ht_l <- lapply(seq_along(EH_l3), function(ihmrow){
         HM_temp <- Reduce("+", EH_l3[[ihmrow]]);
         main_heatmap_temp <- main_heatmap;
         if (length(partition) > 0) {
            HM_temp <- PHM[match(rows, PHM_rows),] + HM_temp;
            main_heatmap_temp <- main_heatmap_temp + 1;
         }
         if (length(AHM) > 0) {
            HM_temp <- AHM[match(rows, AHM_rows),] + HM_temp;
            main_heatmap_temp <- main_heatmap_temp + 1;
         }
         extra_legend <- NULL;
         # add custom legend to the last row
         if (FALSE && ihmrow == length(EH_l3)) {
            extra_legend <- c(
               caption_legends,
               list(top_legend))
            # caption_legendlist);
         }
         ht_1 <- grid::grid.grabExpr(
            ComplexHeatmap::draw(HM_temp,
               ht_gap=ht_gap,
               adjust_annotation_extension=TRUE,
               annotation_legend_list=extra_legend,
               main_heatmap=main_heatmap_temp));
         ht_1;
      });
      if (do_plot) {
         l <- grid::grid.layout(hm_nrow, 1);
         vp <- grid::viewport(width=1, height=1, layout=l);
         grid::grid.newpage();
         grid::pushViewport(vp);
         for (i in seq_along(ht_l)) {
            grid::pushViewport(grid::viewport(layout.pos.row=i));
            grid::grid.draw(ht_l[[i]]);
            grid::popViewport();
         }
         grid::popViewport();
      }
      if ("grid" %in% return_type) {
         EH_l <- ht_l;
      }
   } else {
      ################################
      ## Single row layout
      HM_temp <- Reduce("+", EH_l);
      main_heatmap_temp <- head(main_heatmap, 1);
      ## test to force first heatmap to have row labels
      main_heatmap_temp <- 1;

      ht_gap <- rep(ht_gap,
         length.out=max(c(1, length(nmatlist)-1)));
      if (length(panel_groups) > 0) {
         ht_gap_adjust <- head(
            (panel_groups != tail(c(panel_groups, "blahblah"), -1)) * 2.5 + 0.5,
            -1);
         ht_gap <- ht_gap * ht_gap_adjust;
      }
      if (length(partition) > 0) {
         HM_temp <- PHM[match(rows, PHM_rows),] + HM_temp;
         main_heatmap_temp <- main_heatmap_temp + 1;
         ht_gap <- grid::unit.c(grid::unit(1, "mm"), ht_gap);
      }
      if (length(AHM) > 0) {
         HM_temp <- AHM[match(rows, AHM_rows),] + HM_temp;
         main_heatmap_temp <- main_heatmap_temp + 1;
         ht_gap <- grid::unit.c(grid::unit(1, "mm"), ht_gap);
      }
      if (length(MHM) > 0) {
         HM_temp <- HM_temp + MHM[match(rows, MHM_rows),];
         ht_gap <- grid::unit.c(ht_gap, grid::unit(1, "mm"));
      }
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "ht_gap:");
         print(ht_gap);
      }
      if (do_plot &&
            (length(title) > 0 || length(caption) > 0)) {
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Preparing ComplexHeatmap::draw(HeatmapList)");
         }
         HM_drawn <- ComplexHeatmap::draw(HM_temp,
            column_title=title,
            column_title_gp=title_gp,
            ht_gap=ht_gap,
            adjust_annotation_extension=TRUE,
            annotation_legend_list=c(
               caption_legends,
               list(top_legend)),
            main_heatmap=main_heatmap_temp,
            merge_legends=TRUE,
            padding=padding)
         if (FALSE) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Preparing HeatmapList grob for grid_with_title()");
            HM_grob <- grid::grid.grabExpr(
               ComplexHeatmap::draw(HM_temp,
                  ht_gap=ht_gap,
                  main_heatmap=main_heatmap_temp)
            );
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Calling grid_with_title()");
            multienrichjam::grid_with_title(HM_grob,
               title=title,
               caption=caption,
               verbose=verbose,
               ...);
         }
      } else {
         if (do_plot) {
            HM_drawn <- ComplexHeatmap::draw(HM_temp,
               ht_gap=ht_gap,
               adjust_annotation_extension=TRUE,
               annotation_legend_list=list(top_legend),
               main_heatmap=main_heatmap_temp,
               merge_legends=TRUE,
               padding=padding)
         }
      }
   }

   # revert custom row/column heatmap annotation padding
   ComplexHeatmap::ht_opt("ROW_ANNO_PADDING"=ROW_ANNO_PADDING)
   ComplexHeatmap::ht_opt("COLUMN_ANNO_PADDING"=COLUMN_ANNO_PADDING)
   ComplexHeatmap::ht_opt("HEATMAP_LEGEND_PADDING"=HEATMAP_LEGEND_PADDING)
   ComplexHeatmap::ht_opt("ANNOTATION_LEGEND_PADDING"=ANNOTATION_LEGEND_PADDING)

}
