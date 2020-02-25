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

   ## Interpret vector of filenames as input
   if (length(filename) == 0 && length(x) > 1 && is.atomic(x)) {
      if (all(file.exists(x))) {
         xseq <- nameVector(seq_along(x), names(x));
         nmatlist <- lapply(xseq, function(i){
            coverage_matrix2nmat(x=x[i],
               filename=NULL,
               signal_name=names(x)[i],
               target_name=target_name,
               background=background,
               smooth=smooth,
               target_is_single_point=target_is_single_point,
               signal_is_categorical=signal_is_categorical,
               mat_grep=mat_grep,
               upstream_grep=upstream_grep,
               downstream_grep=downstream_grep,
               target_grep=target_grep,
               verbose=verbose,
               ...);
         });
         return(nmatlist);
      } else {
         stop("When x is supplied as a vector, it should contain files accessible with file.exists(x).");
      }
   }

   ## Other single-entry input types
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
   if (verbose) {
      jamba::printDebug("matrix2nmat(): ",
         "signal_name:",
         signal_name);
      jamba::printDebug("matrix2nmat(): ",
         "names(filename):",
         names(filename));
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
#' @param row_order integer vector used to order rows.
#'    When `TRUE` or `NULL` it uses
#'    the default for `EnrichedHeatmap::EnrichedHeatmap()`
#'    which is the `EnrichedHeatmap::enriched_score()`
#'    for the matrix `main_heatmap`. When `FALSE` the
#'    rows are ordered by the order they appear in `rows`,
#'    which is either the order they appear in `nmatlist`
#'    or the order after sorting `anno_df`. When
#'    `TRUE` the default
#' @param nmat_colors named character vector of R colors,
#'    to colorize each heatmap. When `NULL` then
#'    `colorjam::rainbowJam()` is used to create colors
#'    for each heatmap panel.
#' @param main_heatmap integer index referring to the
#'    entry in `nmatlist` to use for clustering and row
#'    ordering.
#' @param anno_df `data.frame` or object that can be coerced,
#'    used to annotate rows of each matrix. It must have
#'    `rownames(anno_df)` that match `rownames(nmatlist)`.
#'    When supplied, data can be sorted using `byCols`.
#'    Note that only the `rownames(anno_df)`
#'    present in both `nmatlist` and `anno_df` are
#'    used to display the heatmaps. These rows
#'    may also be subsetted using argument `rows`.
#' @param byCols character vector of  values in
#'    `colnames(anno_df)` used to sort the data.frame
#'    via `jamba::mixedSortDF()`. Any colname with
#'    prefix `-` will be reverse-sorted.
#' @param anno_row_marks character vector of `rownames`
#'    which will be labeled beside the heatmaps, using
#'    the `ComplexHeatmap::anno_mark()` method. It currently
#'    requires `anno_df` be defined, since it uses the
#'    first column in `anno_df` as a one-column heatmap,
#'    to anchor the labels.
#' @param anno_row_labels character vector of optional
#'    character labels to use instead of `rownames`.
#'    If `NULL` then `anno_row_marks` are used. Or
#'    `anno_row_labels` may contain a character vector
#'    of `colnames(anno_df)` which will create labels
#'    by concatenating each column value separated by
#'    space `" "`.
#' @param hm_nrow integer number of rows used to display
#'    the heatmap panels.
#' @param transform either `character` string referring to
#'    a numeric transformation, or a `function` that applies
#'    a numeric transformation. Valid `character` string values:
#'    `"log2signed"` applies `jamba::log2signed()` which applies
#'    `log2(1+x)` transform to the absolute value, then multiplies
#'    by the original `sign(x)`; `"sqrt"` applies square root;
#'    `"cubert"` applies cube root `x^(1/3)`; `"qrt"` applies
#'    fourth root `x^(1/4)`. When there are negative numeric
#'    values, the transformation is applied to absolute value,
#'    then multiplied by the original sign. Therefore, the
#'    transformation is applied to adjust the magnitude of
#'    the values. These values are passed to `get_numeric_transform()`
#'    which may have more information.
#' @param signal_ceiling numeric vector length `length(nmatlist)`
#'    which applies a maximum numeric value to the
#'    color ramp for each matrix in `nmatlist`. Every numeric
#'    value above the `signal_ceiling` will be assigned the
#'    maximum color. The values in `signal_ceiling` are
#'    recycled to `length(nmatlist)`, so if one value is
#'    provided, it will be applied to every matrix.
#' @param lens numeric value used to scale each heatmap
#'    color ramp, using `getColorRamp()`. Values above zero
#'    apply the color gradient more rapidly starting from the
#'    lowest value, making the color appear more intense for
#'    lower numeric values. Values below zero apply the color gradient
#'    less rapidly, which makes lower numeric values appear
#'    less intense. This adjustment is intended to help
#'    apply suitable color contrast depending upon the range
#'    of numeric values. The `lens` values are applied to
#'    each matrix in `nmatlist`, and so it is recycled to
#'    `length(nmatlist)` as needed. Note that `signal_ceiling`
#'    is also intended to help apply the color gradient to
#'    a suitable numeric range, and the `lens` argument is
#'    applied relative to the numeric range being used.
#' @param axis_name_gp x-axis label graphic parameters,
#'    as output from `grid::gpar()`. For example to define
#'    the x-axis font size, use the form
#'    `grid::gpar(fontsize=8)`.
#' @param  axis_name_rot numeric value either `0` or `90` indicating
#'    whether to rotate the x-axis names, where `90` will rotate
#'    labels, and `0` will leave labels horizontal.
#' @param seed numeric value used with `set.seed()` to
#'    set the random seed. Set to `NULL` to avoid running
#'    `set.seed()`.
#' @param use_raster logical indicating whether to create heatmaps
#'    using raster resizing, almost always recommended `TRUE`.
#' @param do_plot logical indicating whether to draw the heatmaps,
#'    where `FALSE` will return the data used to create heatmaps
#'    without actually drawing the heatmaps.
#' @param return_type character string indicating the type of
#'    data to return: `"heatmaplist"` returns the list of heatmaps,
#'    which can separately be arranged together using
#'    `ComplexHeatmap::draw()` or `grid::grid.draw()`.
#' @param show_error logical indicating whether to add error
#'    bars to the profile plot at the top of each heatmap.
#'    These error bars are calculated by
#'    `EnrichedHeatmap::anno_enriched()` using
#'    `matrixStats::colSds(x)/nrow(x)`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are passed to
#'    `EnrichedHeatmap::EnrichedHeatmap()` to allow greater
#'    customization of details. Note that many `...` arguments
#'    are also passed to `ComplexHeatmap::Heatmap()`.
#'
#' @family jam genome functions
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
 k_method=c("euclidean", "pearson", "correlation"),
 partition=NULL,
 rows=NULL,
 row_order=NULL,
 nmat_colors=NULL,
 main_heatmap=1,
 anno_df=NULL,
 byCols=NULL,
 anno_row_marks=NULL,
 anno_row_labels=NULL,
 hm_nrow=1,
 transform=jamba::log2signed,
 signal_ceiling=NULL,
 axis_name_gp=grid::gpar(fontsize=8),
 axis_name_rot=90,
 lens=-2,
 seed=123,
 use_raster=TRUE,
 do_plot=TRUE,
 return_type=c("heatmaplist", "grid"),
 show_error=FALSE,
 verbose=FALSE,
 ...)
{
   #
   return_type <- match.arg(return_type);
   if (length(seed) > 0) {
      set.seed(seed);
   }
   if (length(main_heatmap) == 0 || main_heatmap > length(nmatlist)) {
      main_heatmap <- 1;
   }
   ## k_method
   kmeans <- stats::kmeans;
   k_method <- head(k_method, 1);
   if (jamba::igrepHas("pearson|correlation|spearman", k_method)) {
      if (!suppressPackageStartupMessages(require(amap))) {
         k_method <- "euclidean";
         jamba::printDebug("nmatlist2heatmaps(): ",
            "k_method requires the ",
            "amap",
            " package, which is not installed. Setting k_method to ",
            "euclidean");
      }
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "Using amap::Kmeans()");
      }
      kmeans <- function(...){amap::Kmeans(..., method=k_method)};
   }
   if (length(k_method) == 0 || nchar(k_method) == 0) {
      k_method <- "euclidean";
   }
   if (length(rows) == 0) {
      ## Make sure rows are present in all nmatlist entries.
      rows <- Reduce("intersect",
         lapply(nmatlist, rownames));
      #rows <- rownames(nmatlist[[main_heatmap]]);
   } else {
      ## Make sure rows are present in all rownames of nmatlist
      rows <- Reduce("intersect",
         c(list(rows), lapply(nmatlist, rownames)));
   }
   if (length(nmat_colors) == 0) {
      nmat_colors <- colorjam::rainbowJam(length(nmatlist),
         ...);
   }
   if (length(nmat_colors) < length(nmatlist)) {
      nmat_colors <- rep(nmat_colors,
         length.out=length(nmatlist));
   }

   ## Optional transformation of each matrix
   if (length(transform) == 0) {
      transform <- function(x){x}
   }
   transform <- get_numeric_transform(transform);
   if (!is.list(transform)) {
      transform <- list(transform);
   }
   if (length(transform) != length(nmatlist)) {
      transform <- rep(transform,
         length.out=length(nmatlist));
   }
   if (verbose) {
      jamba::printDebug("nmatlist2heatmaps(): ",
         "str(transform):");
      print(str(transform));
   }

   ## optional signal_ceiling
   if (length(signal_ceiling) > 0) {
      signal_ceiling <- rep(signal_ceiling,
         length.out=length(nmatlist));
   }

   ## Define some empty variables
   PHM <- NULL;
   AHM <- NULL;
   MHM <- NULL;

   ## Optional data.frame with additional annotations
   if (length(anno_df) > 0) {
      if (!igrepHas("data.frame|dataframe|data.table|tibble", class(anno_df))) {
         stop("anno_df must be data.frame, DataFrame, data.table, or tibble.");
      }
      if (!any(rownames(anno_df) %in% rows)) {
         stop("anno_df must contain rownames present in nmatlist.");
      }
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "Preparing anno_df.");
      }
      if (length(byCols) > 0) {
         anno_df <- jamba::mixedSortDF(anno_df, byCols=byCols);
         rows <- intersect(rownames(anno_df), rows);
         row_order <- FALSE;
      } else {
         rows <- intersect(rows, rownames(anno_df));
      }
      anno_df <- anno_df[rows,,drop=FALSE];
      ## Determine a list of color functions, one for each column
      anno_colors_l <- lapply(nameVector(colnames(anno_df)), function(i){
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "anno_colors_l colname:", i);
         }
         i1 <- anno_df[[i]];
         if (any(c("integer", "numeric") %in% class(i1))) {
            if (min(i1, na.rm=TRUE) < 0) {
               ## Bi-directional color scale
               #ibreaks1 <- max(abs(i1), na.rm=TRUE);
               ibreaks1 <- quantile(abs(i1), c(0.995));
               ibreaks <- unique(seq(from=-ibreaks1,
                  to=ibreaks1,
                  length.out=25));
               if (length(ibreaks) <= 1) {
                  ibreaks <- c(-1, 0, 1);
               }
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "anno_colors_l colname:", i,
                     " bi-directional data");
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "ibreaks:", ibreaks);
               }
               cBR <- circlize::colorRamp2(breaks=ibreaks,
                  col=jamba::getColorRamp("RdBu_r",
                     lens=2,
                     n=length(ibreaks)));
            } else {
               #ibreaks1 <- max(abs(i1), na.rm=TRUE);
               #ibreaks2 <- min(abs(i1), na.rm=TRUE);
               #imin <- max(c(0, ibreaks2-(ibreaks1 - ibreaks2)*0.2));
               #ibreaks <- seq(from=imin, to=ibreaks1, length.out=15);
               iminmax <- quantile(rmNA(i1),
                  c(0.005, 0.995));
               ibreaks <- unique(seq(from=iminmax[1],
                  to=iminmax[2],
                  length.out=15));
               if (max(ibreaks) == 0 || length(ibreaks) <= 1) {
                  ibreaks <- c(0, 1);
               }
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "anno_colors_l colname:", i,
                     " bi-directional data");
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "ibreaks:", ibreaks);
               }
               cBR <- circlize::colorRamp2(breaks=ibreaks,
                  col=jamba::getColorRamp("Purples",
                     n=length(ibreaks),
                     lens=2));
            }
         } else {
            i2 <- mixedSort(unique(i1));
            if (!"factor" %in% class(i2)) {
               i2 <- factor(i2, levels=mixedSort(i2));
            }
            cBR <- colorjam::group2colors(i2);
         }
         cBR;
      });
      for (jj in 1:ncol(anno_df)) {
         i1 <- anno_df[[jj]];
         if (any(c("character") %in% class(i1))) {
            i1 <- factor(i1, levels=mixedSort(unique(i1)));
            anno_df[[jj]] <- i1;
         }
      }
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "Creating AHM.");
      }
      AHM <- rowAnnotation(df=anno_df[rows,,drop=FALSE],
         name="Annotation",
         col=anno_colors_l);

      ## Optional row marks
      anno_rows <- intersect(rows, anno_row_marks);
      if (length(anno_rows) > 0) {
         anno_row_which <- match(anno_rows, rows);
         if (length(anno_row_labels) > 0 && all(anno_row_labels %in% colnames(anno_df))) {
            anno_row_labels <- pasteByRow(
               anno_df[anno_rows,anno_row_labels,drop=FALSE],
               sep=" ");
         } else if (length(anno_row_labels) >= length(anno_rows)) {
            anno_row_labels <- anno_row_labels[anno_rows];
         } else {
            anno_row_labels <- anno_rows;
         }
         ## Print optional verbose output
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Preparing row mark MHM, top 20 entries are shown:");
            print(head(
               data.frame(
                  anno_rows=anno_rows,
                  anno_row_which=anno_row_which,
                  anno_row_labels=anno_row_labels),
               20));
         }
         ## Mark Heatmap
         MHM <- Heatmap(nameVector(anno_df[rows,1], rows),
            split=partition[rows],
            col=anno_colors_l[[1]],
            use_raster=use_raster,
            name=colnames(anno_df)[1],
            show_row_names=FALSE,
            width=k_width,
            cluster_rows=FALSE,
            right_annotation=rowAnnotation(
               foo=anno_mark(at=anno_row_which,
                  labels=anno_row_labels)
            )
         );
      } else {
         MHM <- NULL;
      }
   } else {
      AHM <- NULL;
      MHM <- NULL;
   }
   if (length(row_order) == 0 | isTRUE(row_order)) {
      row_order <- order(enriched_score(nmatlist[[main_heatmap]][rows,]),
         decreasing=TRUE);
   } else if (isFALSE(row_order)) {
      row_order <- seq_along(rows);
   }

   ## Optional k-means clustering
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
      itransform <- transform[[main_heatmap]];
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "Running kmeans, k_clusters:",
            k_clusters,
            ", k_method:",
            k_method);
      }
      partition <- kmeans(
         itransform(nmatlist[[main_heatmap]][rows,]),
         iter.max=20,
         centers=k_clusters)$cluster;
      if (verbose) {
         jamba::printDebug("table(partition):");
         print(table(partition));
      }
   }
   ## Partition heatmap sidebar
   if (length(partition) > 0) {
      ## Define colors if not provided
      if (length(k_colors) == 0) {
         k_colors <- nameVector(
            colorjam::rainbowJam(length(unique(partition)),
               ...),
            unique(partition));
         k_colors <- k_colors[sort(names(k_colors))];
      }

      ## Optional subset of k-means clusters
      if (length(k_subset) > 0) {
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Subsetting partition for k_subset:",
               k_subset);
         }
         partition <- partition[partition %in% k_subset];
         rows <- names(partition);
         k_colors <- k_colors[names(k_colors) %in% as.character(k_subset)];
         ## Subset AHM and MHM if defined
         if (length(AHM) > 0) {
            AHM <- AHM[match(rows, rownames(attr(AHM, "matrix"))),];
         }
         if (length(MHM) > 0) {
            MHM <- MHM[match(rows, rownames(attr(MHM, "matrix"))),];
         }
      }
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "creating partition heatmap PHM.");
      }
      PHM <- Heatmap(partition[rows],
         split=partition[rows],
         use_raster=use_raster,
         col=k_colors,
         name="cluster",
         show_row_names=FALSE,
         width=k_width);
   }

   ## Iterate each matrix to create heatmaps
   lens <- rep(lens, length.out=length(nmatlist));
   if ("gpar" %in% class(axis_name_gp)) {
      axis_name_gp <- rep(list(axis_name_gp),
         length.out=length(nmatlist))
   } else if (!is.list(axis_name_gp) && length(axis_name_gp) > 0) {
      axis_name_gp <- rep(list(axis_name_gp), length.out=length(nmatlist));
   } else {
      axis_name_gp <- rep(axis_name_gp, length.out=length(nmatlist));
   }
   EH_l <- lapply(seq_along(nmatlist), function(i){
      nmat <- nmatlist[[i]][rows,,drop=FALSE];
      signal_name <- attr(nmat, "signal_name");
      target_name <- attr(nmat, "target_name");
      s_name <- gsub("_at_", "\nat_", signal_name);
      color <- nmat_colors[[i]];
      if (length(color) == 0 || is.na(color)) {
         color <- "aquamarine4";
      }
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "signal_name:",
            signal_name,
            ", target_name:",
            target_name,
            ", color:",
            color,
            fgText=c(
               rep(
                  list("darkorange",
                     "dodgerblue"),
                  length.out=6),
               color));
      }
      itransform <- transform[[i]];
      imat <- itransform(nmat);
      iceiling <- signal_ceiling[[i]];
      if (length(iceiling) > 0 && !is.na(iceiling)) {
         if (iceiling > 0 && iceiling <= 1) {
            # apply quantile
            imat_values <- setdiff(abs(imat), 0);
            iquantile <- quantile(imat_values,
               probs=iceiling,
               na.rm=TRUE);
            if (verbose) {
               jamba::printDebug("nmatlist2heatmaps(): ",
                  "applied iceiling:",
                  iceiling,
                  " as quantile threshold, which defined new ceiling:",
                  iquantile);
            }
            iceiling <- iquantile;
         }
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "   Applied ceiling:",
               iceiling);
         }
         if (min(imat, na.rm=TRUE) < 0) {
            ibreaks <- seq(from=-iceiling,
               to=iceiling,
               length=21);
         } else {
            ibreaks <- seq(from=0,
               to=iceiling,
               length=21);
         }
         colramp <- circlize::colorRamp2(
            breaks=ibreaks,
            colors=jamba::getColorRamp(color,
               n=21,
               lens=lens[[i]]));
      } else {
         colramp <- jamba::getColorRamp(color,
            n=21,
            lens=lens[[i]]);
      }
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "axis_name_gp[[i]]:");
         print(axis_name_gp[[i]]);
      }
      EH <- EnrichedHeatmap::EnrichedHeatmap(imat,
         split=partition[rows],
         pos_line=FALSE,
         use_raster=use_raster,
         col=colramp,
         top_annotation=ComplexHeatmap::HeatmapAnnotation(
            lines=EnrichedHeatmap::anno_enriched(gp=grid::gpar(col=k_colors),
               show_error=show_error)
         ),
         axis_name_gp=axis_name_gp[[i]],
         axis_name_rot=axis_name_rot,
         name=signal_name,
         column_title=signal_name,
         row_order=row_order,
         ...);
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
         HM_temp <- Reduce("+", EHs);
         main_heatmap_temp <- main_heatmap;
         if (length(partition) > 0) {
            HM_temp <- PHM + HM_temp;
            main_heatmap_temp <- main_heatmap_temp + 1;
         }
         if (length(AHM) > 0) {
            HM_temp <- AHM + HM_temp;
            main_heatmap_temp <- main_heatmap_temp + 1;
         }
         ht_1 <- grid::grid.grabExpr(
            ComplexHeatmap::draw(HM_temp,
               main_heatmap=main_heatmap_temp));
         ht_1;
      });
      if (do_plot) {
         l <- grid::grid.layout(hm_nrow, 1);
         vp <- grid::viewport(width=1, height=1, layout=l);
         grid::grid.newpage();
         grid::pushViewport(vp);
         for (i in seq_along(ht_l)) {
            grid::pushViewport(viewport(layout.pos.row=i));
            grid::grid.draw(ht_l[[i]]);
            grid::popViewport();
         }
         grid::popViewport();
      }
      if ("grid" %in% return_type) {
         EH_l <- ht_l;
      }
   } else {
      ## Single row layout
      HM_temp <- Reduce("+", EH_l);
      main_heatmap_temp <- main_heatmap;
      if (length(partition) > 0) {
         HM_temp <- PHM + HM_temp;
         main_heatmap_temp <- main_heatmap_temp + 1;
      }
      if (length(AHM) > 0) {
         HM_temp <- AHM + HM_temp;
         main_heatmap_temp <- main_heatmap_temp + 1;
      }
      if (length(MHM) > 0) {
         HM_temp <- HM_temp + MHM;
      }
      ComplexHeatmap::draw(HM_temp,
         main_heatmap=main_heatmap_temp);
   }
   invisible(c(list(AHM=AHM),
      list(PHM=PHM),
      EH=EH_l,
      list(MHM=MHM)));
}


#' Import deepTools coverage matrix to normalizedMatrix
#'
#' Import deepTools coverage matrix to normalizedMatrix
#'
#' This function is under active development.
#'
#' @export
deepTools_matrix2nmat <- function
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
   x1 <- readLines(filename, n=1);
   xyaml <- yaml::read_yaml(text=gsub("^@", "", x1));
   samples <- xyaml$sample_labels;
   samples;
   binsize <- xyaml$`bin size`;
   binsize;
   upstream <- xyaml$upstream;
   upstream;
   group_boundaries <- xyaml$group_boundaries;
   sample_boundaries <- xyaml$sample_boundaries;
   x <- data.table::fread(text=readLines(filename),
      skip=1,
      header=FALSE,
      sep="\t");
   rownames(x) <- x[[4]];

   mats <- as.matrix(x[,-1:-6,drop=FALSE]);
   mat_n <- ceiling(ncol(mats) / length(samples));

   ## required attributes
   starts <- seq(from=-upstream, by=binsize, length.out=mat_n);
   upstream_index <- which(starts < 0);
   target_index <- which(starts == 0);
   downstream_index <- which(starts > 0);
   target_name <- xyaml$group_labels;

   mat_split <- rep(rep(samples, each=mat_n), length.out=ncol(mats));
   colnames_l <- split(colnames(mats), mat_split);
   mat_l <- lapply(seq_along(colnames_l), function(k){
      i <- colnames_l[[k]];
      im <- mats[,i];
      colnames(im) <- paste0("u", seq_len(ncol(im)));
      rownames(im) <- x[[4]];
      attr(im, "signal_name") <- samples[k];
      attr(im, "target_name") <- target_name;
      attr(im, "upstream_index") <- upstream_index;
      attr(im, "target_index") <- target_index;
      attr(im, "downstream_index") <- downstream_index;
      im;
   });

}


#' Get appropriate numeric transformation function
#'
#' Get appropriate numeric transformation function
#'
#' This function recognizes numeric transformation functions by
#' name, or searches the attached R package environments for
#' a matching function by name.
#'
#' @return `function` or `NULL` when no matching function is
#'    found, or `list` is returned when the input `transform`
#'    has multiple values.
#'
#' @param transform `character` string or `function`. If
#'    `length(transform) > 1` then `get_numeric_transform()`
#'    is called for each value in `transform`. Valid names:
#'    `"log2signed"` applies `jamba::log2signed()` which applies
#'    `log2(1+x)` transform to the absolute value, then multiplies
#'    by the original `sign(x)`;
#'    `"exp2signed"` applies the inverse of `"log2signed"`, which
#'    exponentiates numeric values which were previously transformed
#'    with `log2(1+x)`;
#'    `"sqrt"` applies square root transform;
#'    `"cubert"` applies cube root `x^(1/3)`;
#'    `"qrt"` applies fourth root `x^(1/4)`;
#'    `"frt"` applies fourth root `x^(1/5)`;
#'    `"square"` applies `x^2` to absolute value, multiplied by the `sign(x)`;
#'    `"cube"` applies `x^3`;
#'    When there are negative numeric
#'    values, the transformation is applied to absolute value,
#'    then multiplied by the original sign. Therefore, the
#'    transformation is applied to adjust the magnitude of
#'    the values.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' get_numeric_transform("log2signed")
#'
#' @export
get_numeric_transform <- function
(transform,
 ...)
{
   #
   if (length(transform) > 1) {
      transform <- lapply(transform, function(it){
         get_numeric_transform(it, ...)
      });
      return(transform);
   }
   it <- transform;
   if (is.atomic(it)) {
      if ("log2signed" %in% it) {
         it <- jamba::log2signed;
      } else if ("exp2signed" %in% it) {
         it <- jamba::exp2signed;
      } else if ("sqrt" %in% it) {
         it <- function(x){sign(x)*sqrt(x)};
      } else if ("square" %in% it) {
         it <- function(x){sign(x)*(abs(x)^2)};
      } else if ("cubert" %in% it) {
         it <- function(x){sign(x)*abs(x)^(1/3)};
      } else if ("cube" %in% it) {
         it <- function(x){sign(x)*(abs(x)^3)};
      } else if ("qrt" %in% it) {
         it <- function(x){sign(x)*abs(x)^(1/4)};
      } else if ("frt" %in% it) {
         it <- function(x){sign(x)*abs(x)^(1/5)};
      } else {
         itf <- tryCatch({
            get(it);
         }, error=function(e){
            jamba::printDebug("Error:",
               "transform name '",
               it,
               "' was not recognized, and not available on the search path:\n",
               search(),
               sep=", ")
            NULL;
         });
         if (is.function(itf)) {
            it <- function(x){sign(x)*itf(f)};
         } else {
            it <- NULL;
         }
      }
   }
   if (!is.function(it)) {
      return(NULL);
   }
   return(it);
}
