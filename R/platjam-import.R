#' Import genome coverage matrix files
#'
#' Import genome coverage matrix files
#'
#' This function imports genome coverage data matrix
#' and returns an object of class
#' `normalizedMatrix` compatible for use by the
#' package `"EnrichedHeatmap"`.
#'
#' There is a conversion function `EnrichedHeatmap::as.normalizedMatrix()`,
#' however this function does not call that function, in
#' favor of defining the attributes directly. In future, this
#' function may change to call that function.
#'
#' @family jam genome functions
#' @family jam import functions
#'
#' @return `normalizedMatrix` numeric matrix, where additiona
#'    metadata is stored in the object attributes. See
#'    `EnrichedHeatmap::as.normalizedMatrix()` for more
#'    details about the metadata. The `rownames` are defined
#'    by the first colname which does not match
#'    `mat_grep`, which by default is `"Gene ID"`,
#'    otherwise rownames are `NULL`.
#'
#' @param x `data.frame` or compatible object containing
#'    genome coverage data, or a character file path. When
#'    `x` is not supplied, `filename` is used to import
#'    data. When `x` is a filename, it is used to populate
#'    `filename`, then data is imported into `x`.
#' @param filename character path to a genome coverage file.
#'    When `x` is supplied, this argument is ignored. When
#'    `filename` is used, only the first file is imported.
#' @param signal_name The name of signal regions. It is only used
#'    for printing the object. When `signal_name` is `NULL`, the
#'    `signal_name` is derived from `names(filename)` if
#'    available, then `basename(filename)`, or `"signal"` then
#'    only `x` is supplied.
#' @param target_name The name of the target names. It is only used
#'    for printing the object.
#' @param background numeric value containing the background
#'    value in the matrix.
#' @param smooth logical whether to apply smoothing on rows.
#' @param target_is_single_point,signal_is_categorical logical
#'    indicating whether the target region is a single point,
#'    and whether signal matrix is categorical, respectively.
#' @param mat_grep character regular expression pattern used
#'    to identify colnames which contain coverage data. The
#'    default pattern expects the format `"-200:-100"`.
#' @param upstream_grep character regular expression pattern
#'    used to identify upstream colnames from values that
#'    match `mat_grep`. The default assumes any region
#'    beginning `"-"` is negative and upstream the central
#'    target region.
#' @param downstream_grep character regular expression pattern
#'    used to identify upstream colnames from values that
#'    match `mat_grep`. The default assumes all colnames which
#'    are not upstream are therefore downstream.
#' @param target_grep character regular expression pattern
#'    used to identify a colname referring to the `target`,
#'    which by default can only be `"0"`. Otherwise, no target
#'    region is defined.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' ## There is a small example file to use for testing
#' cov_file <- system.file("data", "tss_coverage.matrix", package="platjam");
#' cov_file <- system.file("data", "h3k4me1_coverage.matrix", package="platjam");
#' if (length(cov_file) > 0) {
#'    nmat <- coverage_matrix2nmat(cov_file);
#'    jamba::printDebug("signal_name: ",
#'       attr(nmat, "signal_name"));
#'
#' if (suppressPackageStartupMessages(require(EnrichedHeatmap))) {
#'    color <- "red3";
#'    signal_name <- attr(nmat, "signal_name");
#'    k <- 6;
#'    set.seed(123);
#'    partition <- kmeans(log10(1+nmat), centers=k)$cluster;
#'    EH <- EnrichedHeatmap(log10(1+nmat),
#'       split=partition,
#'       pos_line=FALSE,
#'       use_raster=TRUE,
#'       col=jamba::getColorRamp(color, n=10),
#'       top_annotation=HeatmapAnnotation(
#'          lines=anno_enriched(gp=gpar(col=colorjam::rainbowJam(k)))
#'       ),
#'       axis_name_gp=gpar(fontsize=8),
#'       name=signal_name,
#'       column_title=signal_name
#'    );
#'    PHM <- Heatmap(partition,
#'       use_raster=TRUE,
#'       col=structure(colorjam::rainbowJam(k),
#'          names=as.character(seq_len(k))),
#'       name="partition",
#'       show_row_names=FALSE,
#'       width=unit(3, "mm"));
#'    draw(PHM + EH, main_heatmap=2);
#' }
#' }
#'
#' @export
coverage_matrix2nmat <- function
(x=NULL,
 filename=NULL,
 signal_name=NULL,
 target_name="target",
 background=0,
 smooth=FALSE,
 target_is_single_point=FALSE,
 signal_is_categorical=FALSE,
 mat_grep="[-0-9]+:[-0-9]+",
 upstream_grep="^[-]",
 downstream_grep="^[^-]",
 target_grep="^0$",
 verbose=FALSE,
 ...)
{
   ## read a matrix file, and convert to format EnrichedHeatmap understands
   if (length(x) == 0) {
      if (length(filename) == 0) {
         stop("Must suppled x as data.frame, or filename.");
      } else if (!file.exists(head(x, 1))) {
         stop("No x supplied as data.frame, and filename is not found.");
      } else {
         if (verbose) {
            jamba::printDebug("matrix2nmat(): ",
               "Importing data from filename:",
               filename);
         }
         x <- as.data.frame(
            data.table::fread(head(filename, 1),
               sep="\t"));
         if (length(signal_name) == 0) {
            if (length(names(filename)) > 0 && nchar(head(names(filename), 1)) > 0) {
               signal_name <- head(names(filename), 1);
            } else {
               signal_name <- basename(head(filename, 1));
            }
         }
      }
   } else if (jamba::igrepHas("character", class(x)) && file.exists(head(x, 1))) {
      if (verbose) {
         jamba::printDebug("matrix2nmat(): ",
            "Importing data from filename sent as x:",
            x);
      }
      filename <- x;
      x <- as.data.frame(
         data.table::fread(head(filename, 1),
            sep="\t"));
      if (length(signal_name) == 0) {
         if (length(names(filename)) > 0 && nchar(head(names(filename), 1)) > 0) {
            signal_name <- head(names(filename), 1);
         } else {
            signal_name <- basename(head(filename, 1));
         }
      }
   }

   if (length(signal_name) == 0) {
      signal_name <- "signal";
   }

   if (!jamba::igrepHas("data.frame|dataframe|tibble|data.table", class(x))) {
      stop("Supplied x must be class data.frame, DataFrame, data.table, or tibble.");
   }
   if (jamba::igrepHas("data.table", class(x))) {
      x <- as.data.frame(x);
   }

   ## Extract numeric matrix
   mat_colnames <- jamba::vigrep(mat_grep, colnames(x));
   if (length(mat_colnames) == 0) {
      stop("No coordinate colnames found matching pattern.");
   }
   if (verbose) {
      jamba::printDebug("matrix2nmat(): ",
         "mat_colnames:",
         mat_colnames);
   }
   mat <- as.matrix(x[,mat_colnames,drop=FALSE]);
   name_colname <- head(setdiff(colnames(x), mat_colnames), 1);
   if (length(name_colname) > 0) {
      rownames(mat) <- x[,1];
   }

   ## Extract upstream and downstream
   colnames_up <- jamba::igrep(upstream_grep, mat_colnames);
   colnames_dn <- jamba::igrep(downstream_grep, mat_colnames);
   colnames_target <- jamba::igrep(target_grep, mat_colnames);
   colnames_dn <- setdiff(colnames_dn, colnames_target);
   attr(mat, "upstream_index") <- colnames_up;
   attr(mat, "target_index") <- colnames_target;
   attr(mat, "downstream_index") <- colnames_dn;

   ## Store filename if supplied
   if (length(filename) > 0) {
      attr(mat, "filename") <- filename;
   }

   ## Extract number of bases extended
   mat_extend <- abs(range(as.numeric(unlist(strsplit(mat_colnames, ":")))));
   attr(mat, "extend") <- mat_extend;

   ## Smooth
   attr(mat, "smooth") <- smooth;

   ## Signal and target name
   attr(mat, "signal_name") <- signal_name;
   attr(mat, "target_name") <- target_name;
   attr(mat, "target_is_single_point") <- target_is_single_point;
   attr(mat, "background") <- background;
   attr(mat, "signal_is_categorical") <- signal_is_categorical;
   class(mat) <- c("normalizedMatrix", "matrix");
   mat;
}


#' Make multiple coverage heatmaps
#'
#' Make multiple coverage heatmaps
#'
#' This function takes a list of `normalizedMatrix` objects,
#' usually the output of `coverage_matrix2nmat()`, and
#' produces multiple heatmaps using
#' `EnrichedHeatmap`.
#'
#' This function is intended to be a convenient wrapper to
#' help keep each data matrix in order, to apply consistent
#' clustering and filtering across all data matrices,
#' and to enable optional multi-row heatmap layout.
#'
#' @param nmatlist `list` containing `normalizedMatrix` objects,
#'    usually the output from `coverage_matrix2nmat()`.
#' @param k_clusters integer number of k-means clusters to
#'    use to partition each heatmap. Use `0` or `NULL` for
#'    no clustering.
#' @param k_subset integer vector of k-means clusters to
#'    retain. Often one cluster contains mostly empty
#'    values, and can be removed using this mechanism.
#' @param k_colors vector of R colors, or `NULL` to use
#'    the output of `colorjam::rainbowJam(k_clusters)`.
#' @param k_width unit width of the k-means cluster color
#'    bar, used with `k_clusters`.
#' @param partition vector used to split rows of each
#'    matrix in `nmatlist`, named by rownames. This
#'    value is ignored when `k_clusters` is supplied.
#' @param rows vector of `rownames` or integer vector with
#'    index of rows to keep from each matrix in `nmatlist`.
#' @param nmat_colors named character vector of R colors,
#'    to colorize each heatmap. When `NULL` then
#'    `colorjam::rainbowJam()` is used to create colors
#'    for each heatmap panel.
#' @param main_heatmap integer index referring to the
#'    entry in `nmatlist` to use for clustering and row
#'    ordering.
#' @param hm_nrow integer number of rows used to display
#'    the heatmap panels.
#' @param transform `function` used to transform numeric
#'    values in each entry in `nmatlist`.
#' @param lens numeric value used to scale each heatmap
#'    color ramp, using `getColorRamp()`.
#' @param seed numeric value used with `set.seed()` to
#'    set the random seed. Set to `NULL` to avoid running
#'    `set.seed()`.
#' @param use_raster logical indicating whether to create heatmaps
#'    using raster resizing, almost always recommended `TRUE`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are sent to `colorjam::rainbowJam()`.
#'
#' @examples
#' ## There is a small example file to use for testing
#' cov_file1 <- system.file("data", "tss_coverage.matrix", package="platjam");
#' cov_file2 <- system.file("data", "h3k4me1_coverage.matrix", package="platjam");
#' cov_files <- c(cov_file1, cov_file2);
#' names(cov_files) <- gsub("[.]matrix",
#'    "",
#'    basename(cov_files));
#' nmatlist <- lapply(cov_files, coverage_matrix2nmat);
#' nmatlist2heatmaps(nmatlist);
#'
#' # k-means clusters
#' nmatlist2heatmaps(nmatlist, k_clusters=4);
#'
#' # multiple rows
#' nmatlist2heatmaps(nmatlist, k_clusters=4, hm_nrow=2);
#'
#' @export
nmatlist2heatmaps <- function
(nmatlist,
 k_clusters=0,
 k_subset=NULL,
 k_colors=NULL,
 k_width=unit(3, "mm"),
 partition=NULL,
 rows=NULL,
 nmat_colors=NULL,
 main_heatmap=1,
 hm_nrow=1,
 transform=jamba::log2signed,
 lens=-2,
 seed=123,
 use_raster=TRUE,
 verbose=FALSE,
 ...)
{
   #
   if (length(seed) > 0) {
      set.seed(seed);
   }
   if (length(main_heatmap) == 0 || main_heatmap > length(nmatlist)) {
      main_heatmap <- 1;
   }
   if (length(rows) == 0) {
      rows <- rownames(nmatlist[[main_heatmap]]);
   }
   if (length(nmat_colors) == 0) {
      nmat_colors <- colorjam::rainbowJam(length(nmatlist),
         ...);
   }
   if (length(nmat_colors) < length(nmatlist)) {
      nmat_colors <- rep(nmat_colors,
         length.out=length(nmatlist));
   }
   if (length(k_clusters) > 0 && k_clusters > 0) {
      if (length(k_colors) == 0) {
         k_colors <- jamba::nameVector(
            colorjam::rainbowJam(k_clusters,
               ...),
            seq_len(k_clusters));
      } else if (length(k_colors) < k_clusters) {
         k_colors <- jamba::nameVector(
            rep(k_colors,
               length.out=k_clusters),
            seq_len(k_clusters));
      } else if (length(names(k_colors)) == 0) {
         names(k_colors) <- seq_len(k_clusters);
      }
      partition <- kmeans(
         transform(nmatlist[[main_heatmap]][rows,]),
         centers=k_clusters)$cluster;
      if (length(k_subset) > 0) {
         partition <- partition[partition %in% k_subset];
         rows <- names(partition);
         k_colors <- k_colors[as.character(k_subset)];
      }
   }
   if (length(partition) > 0) {
      if (length(k_colors) == 0) {
         k_colors <- nameVector(
            colorjam::rainbowJam(length(unique(partition)),
               ...),
            unique(partition));
      }
      PHM <- Heatmap(partition,
         use_raster=use_raster,
         col=k_colors,
         name="cluster",
         show_row_names=FALSE,
         width=k_width);
   }
   lens <- rep(lens, length.out=length(nmatlist));
   EH_l <- lapply(seq_along(nmatlist), function(i){
      nmat <- nmatlist[[i]][rows,,drop=FALSE];
      signal_name <- attr(nmat, "signal_name");
      target_name <- attr(nmat, "target_name");
      s_name <- gsub("_at_", "\nat_", signal_name);
      color <- nmat_colors[[i]];
      EH <- EnrichedHeatmap(log10(1+nmat),
         split=partition,
         pos_line=FALSE,
         use_raster=use_raster,
         col=jamba::getColorRamp(color,
            n=10,
            lens=lens[i]),
         top_annotation=HeatmapAnnotation(
            lines=anno_enriched(gp=gpar(col=k_colors))
         ),
         axis_name_gp=gpar(fontsize=8),
         name=signal_name,
         column_title=signal_name
      );
      EH;
   });

   ## Optional multi-row layout
   if (hm_nrow > 1 && length(nmatlist) > 1) {
      hm_split <- rep(rep(
            seq_len(hm_nrow),
            each=ceiling(length(nmatlist) / hm_nrow)),
         length.out=length(nmatlist));
      EH_l3 <- split(EH_l, hm_split);
      ht_l <- lapply(EH_l3, function(EHs){
         if (length(partition) > 0) {
            ht_1 <- grid.grabExpr(
               draw(PHM + Reduce("+", EHs),
                  main_heatmap=main_heatmap+1));
         } else {
            ht_1 <- grid.grabExpr(
               draw(Reduce("+", EHs),
                  main_heatmap=main_heatmap));
         }
         ht_1;
      });
      l <- grid.layout(hm_nrow, 1);
      vp <- viewport(width=1, height=1, layout=l);
      grid.newpage();
      pushViewport(vp);
      for (i in seq_along(ht_l)) {
         pushViewport(viewport(layout.pos.row=i));
         grid.draw(ht_l[[i]]);
         popViewport();
      }
      popViewport();
   } else {
      ## Single row layout
      if (length(partition) > 0) {
         draw(PHM + Reduce("+", EH_l),
            main_heatmap=main_heatmap+1);
      } else {
         draw(Reduce("+", EH_l),
            main_heatmap=main_heatmap);
      }
   }
   invisible(EH_l);
}
