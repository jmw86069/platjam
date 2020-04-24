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
#'          lines=anno_enriched(gp=grid::gpar(col=colorjam::rainbowJam(k)))
#'       ),
#'       axis_name_gp=grid::gpar(fontsize=8),
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
#' @param middle_color `character` R compatible color used
#'    when creating a divergent color gradient, this color
#'    is used as the middle color. Usually this color should
#'    be either `"white"` or `"black"`.
#' @param nmat_names `character` vector, or `NULL`, optional,
#'    used as custom names for each heatmap in `nmatlist`.
#'    When `nmat_names=NULL` the `signal_name` values are
#'    used from each `nmatlist` matrix.
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
#' @param column_title_gp heatmap title graphic parameters,
#'    as output from `grid::gpar()`. For example to define
#'    the x-axis font size, use the form
#'    `grid::gpar(fontsize=8)`. This argument is passed
#'    directly to `ComplexHeatmap::Heatmap()`.
#' @param seed numeric value used with `set.seed()` to
#'    set the random seed. Set to `NULL` to avoid running
#'    `set.seed()`.
#' @param ht_gap unit size to specify the gap between multiple heatmaps.
#'    This argument is passed to `ComplexHeatmap::draw()`. An example
#'    is `grid::unit(8, "mm")` to specify 8 millimeters.
#' @param profile_value character string to define the type of numeric
#'    profile to display at the top of each heatmap. This argument is
#'    passed to `EnrichedHeatmap::anno_enriched()`. Values: `"mean"` the
#'    mean profile; `"sum"` the sum; `"abs_sum"` sum of absolute values;
#'    `"abs_mean"` the mean of absolute values.
#' @param ylims `vector` of maximum y-axis values for each heatmap profile;
#'    or `list`
#' @param border `logical` indicating whether to draw a border around the
#'    heatmap, which includes all heatmap panels in the event of
#'    splitting by clustering. The `border` can be supplied as a vector,
#'    so the `border` can be applied specifically to each heatmap
#'    if needed.
#' @param iter.max integer value indicating the maximum iterations
#'    performed by k-means clustering, only relevant when `k_clusters`
#'    is non-zero.
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
 panel_groups=NULL,
 title=NULL,
 caption=NULL,
 k_clusters=0,
 k_subset=NULL,
 k_colors=NULL,
 k_width=unit(3, "mm"),
 k_method=c("euclidean", "pearson", "correlation"),
 partition=NULL,
 rows=NULL,
 row_order=NULL,
 nmat_colors=NULL,
 middle_color="white",
 nmat_names=NULL,
 main_heatmap=1,
 anno_df=NULL,
 byCols=NULL,
 anno_row_marks=NULL,
 anno_row_labels=NULL,
 hm_nrow=1,
 transform=jamba::log2signed,
 signal_ceiling=NULL,
 axis_name=NULL,
 axis_name_gp=grid::gpar(fontsize=8),
 axis_name_rot=90,
 column_title_gp=gpar(fontsize=12),
 lens=-2,
 seed=123,
 ht_gap=grid::unit(7, "mm"),
 profile_value=c("mean", "sum", "abs_mean", "abs_sum"),
 ylims=NULL,
 border=TRUE,
 iter.max=20,
 use_raster=TRUE,
 do_plot=TRUE,
 legend_width=grid::unit(3, "cm"),
 heatmap_legend_param=NULL,
 annotation_legend_param=NULL,
 return_type=c("heatmaplist", "grid"),
 show_error=FALSE,
 verbose=TRUE,
 ...)
{
   #
   return_type <- match.arg(return_type);
   profile_value <- match.arg(profile_value);
   if (length(seed) > 0) {
      set.seed(seed);
   }
   if (length(main_heatmap) == 0 || main_heatmap > length(nmatlist)) {
      main_heatmap <- 1;
   }

   if (length(border) == 0) {
      border <- FALSE;
   }
   border <- rep(border, length.out=length(nmatlist));
   if (length(legend_width) == 0) {
      legend_width <- grid::unit(3, "cm");
   }

   ## k_method
   kmeans <- stats::kmeans;
   k_method <- head(k_method, 1);
   if (jamba::igrepHas("pearson|correlation|spearman|maximum|manhattan", k_method)) {
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
      kmeans <- function(centers,...){amap::Kmeans(..., centers=centers,method=k_method)};
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
      if (is.numeric(rows)) {
         rows <- Reduce("intersect",
            c(list(rows), lapply(nmatlist, function(im){seq_len(nrow(im))})));
         rows <- rownames(nmatlist[[1]])[rows];
      } else {
         rows <- Reduce("intersect",
            c(list(rows), lapply(nmatlist, rownames)));
      }
   }
   ## row_order must be named by rows
   if (length(row_order) > 1) {
      if (length(names(row_order)) == 0) {
         names(row_order) <- rows;
      }
   }

   if (verbose) {
      jamba::printDebug("nmatlist2heatmaps(): ",
         "Recognized ",
         jamba::formatInt(length(rows)),
         " rows shared across all matrices.");
   }

   if (length(panel_groups) > 0) {
      panel_groups <- rep(panel_groups,
         length.out=length(nmatlist));
   }
   if (length(nmat_colors) == 0) {
      if (length(panel_groups) > 0) {
         nmat_colors <- colorjam::group2colors(panel_groups,
            ...);
      } else {
         nmat_colors <- colorjam::rainbowJam(length(nmatlist),
            ...);
      }
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
   if (verbose > 1) {
      jamba::printDebug("nmatlist2heatmaps(): ",
         "str(transform):");
      print(str(transform));
   }

   ## optional signal_ceiling
   if (length(signal_ceiling) > 0) {
      signal_ceiling <- rep(signal_ceiling,
         length.out=length(nmatlist));
   }
   if (length(ylims) > 0) {
      if (is.list(ylims)) {
         ylims <- rep(ylims,
            length.out=length(nmatlist));
      } else {
         ylims <- lapply(rep(ylims, length.out=length(nmatlist)), function(ylim){
            range(c(0, 0.001, ylim))
         });
      }
   }
   if (length(ylims) > 0 && verbose) {
      jamba::printDebug("nmatlist2heatmaps(): ",
         "recognized ylims: ",
         paste0("(", jamba::cPaste(ylims), ")"),
         sep="; ");
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
      if (!all(rows %in% rownames(anno_df)) && verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "Using subset of ",
            jamba::formatInt(sum(rows %in% rownames(anno_df))),
            " rows present in rownames(anno_df).");
      }
      if (length(byCols) > 0) {
         anno_df <- jamba::mixedSortDF(anno_df,
            byCols=byCols);
         rows <- intersect(rownames(anno_df), rows);
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Sorted rows by:",
               byCols);
         }
         row_order <- FALSE;
      } else {
         rows <- intersect(rows,
            rownames(anno_df));
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
               colBR <- jamba::getColorRamp("RdBu_r",
                  lens=1,
                  n=length(ibreaks));
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "anno_colors_l colname:", i,
                     " bi-directional data");
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "ibreaks:");
                  jamba::printDebugI(
                     jamba::nameVector(colBR, round(digits=2, ibreaks)),
                     sep=", ");
               }
               cBR <- circlize::colorRamp2(breaks=ibreaks,
                  col=colBR);
            } else {
               iminmax <- quantile(rmNA(i1),
                  c(0.005, 0.995));
               ibreaks <- unique(seq(from=iminmax[1],
                  to=iminmax[2],
                  length.out=15));
               if (max(ibreaks) == 0 || length(ibreaks) <= 1) {
                  ibreaks <- c(0, 1);
               }
               colBR <- jamba::getColorRamp("Purples",
                  n=length(ibreaks),
                  lens=1);
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "anno_colors_l colname:", i,
                     " uni-directional data");
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "ibreaks:");
                  jamba::printDebugI(
                     jamba::nameVector(colBR, round(digits=2, ibreaks)),
                     sep=", ");
               }
               cBR <- circlize::colorRamp2(breaks=ibreaks,
                  col=colBR);
            }
         } else {
            i2 <- jamba::mixedSort(unique(i1));
            if (!"factor" %in% class(i2)) {
               i2 <- factor(i2, levels=jamba::mixedSort(i2));
            }
            cBR <- colorjam::group2colors(i2);
         }
         cBR;
      });
      ## annotation_legend_param
      if (length(annotation_legend_param) == 0) {
         annotation_legend_param <- lapply(nameVector(colnames(anno_df)), function(i){
            i1 <- anno_df[[i]];
            if (any(c("integer", "numeric") %in% class(i1))) {
               list(
                  direction="horizontal",
                  title=i,
                  legend_width=legend_width,
                  title_position="topleft",
                  border="black",
                  grid_width=unit(1, "npc"));
            } else {
               list(
                  title=i,
                  title_position="topleft",
                  border="black",
                  ncol=min(c(length(unique(i1)), 4))
               )
               #   grid_width=unit(1, "npc")
            }
         });
      }
      for (jj in 1:ncol(anno_df)) {
         i1 <- anno_df[[jj]];
         if (any(c("character") %in% class(i1))) {
            i1 <- factor(i1, levels=jamba::mixedSort(unique(i1)));
            anno_df[[jj]] <- i1;
         }
      }
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "Creating anno_df heatmap AHM.");
      }
      ##################################
      ## Annotation heatmap
      AHM <- ComplexHeatmap::rowAnnotation(
         df=anno_df[rows,,drop=FALSE],
         annotation_legend_param=annotation_legend_param,
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
               "Preparing row marks for ",
               jamba::formatInt(length(anno_rows)),
               " anno_rows found in anno_df, top 20 entries are shown:");
            print(head(
               data.frame(
                  anno_rows=anno_rows,
                  anno_row_which=anno_row_which,
                  anno_row_labels=anno_row_labels),
               20));
         }
         ##################################
         ## Mark Heatmap
         MHM <- ComplexHeatmap::Heatmap(nameVector(anno_df[rows,1], rows),
            #split=partition[rows],
            col=anno_colors_l[[1]],
            #use_raster=use_raster,
            name=colnames(anno_df)[1],
            show_row_names=FALSE,
            width=k_width,
            cluster_rows=FALSE,
            right_annotation=ComplexHeatmap::rowAnnotation(
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
         iter.max=iter.max,
         centers=k_clusters)$cluster;
      ## Confirm that names(partition) match rows
      names(partition) <- rows;
      if (verbose) {
         k_sizes <- table(partition);
         jamba::printDebug("nmatlist2heatmaps(): ",
            "k-means cluster sizes: ",
            paste0("cluster", names(k_sizes), "=", k_sizes), sep=", ");
      }
   }
   ## Partition heatmap sidebar
   if (length(partition) > 0) {
      ## Make sure to use the partition values with the properly ordered rows
      if (!all(rows %in% names(partition))) {
         jamba::printDebug("head(partition):");
         print(head(partition));
         jamba::printDebug("head(rows):");
         print(head(rows));
         print(table(all(rows) %in% names(partition)));
         stop("names(partition) must match rownames in nmatlist.");
      }
      partition <- partition[rows];
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
         partition <- partition[as.character(partition) %in% as.character(k_subset)];
         rows <- names(partition);
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Subsetting partition for k_subset:",
               k_subset,
               " from ",
               jamba::formatInt(length(k_colors)),
               " rows down to ",
               jamba::formatInt(length(rows)));
         }
         k_colors <- k_colors[names(k_colors) %in% as.character(k_subset)];
         ## Subset AHM and MHM if defined
         if (length(AHM) > 0) {
            AHM <- AHM[rows,];
            if (verbose) {
               jamba::printDebug("nmatlist2heatmaps(): ",
                  "Taking k_subset rows for annotation heatmap.");
            }
         }
         if (length(MHM) > 0) {
            MHM <- MHM[rows,];
            if (verbose) {
               jamba::printDebug("nmatlist2heatmaps(): ",
                  "Taking k_subset rows for annotation mark heatmap.");
            }
         }
      }
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "creating partition heatmap PHM.");
      }
      ##################################
      ## Partition Heatmap
      p_annotation_legend_param <- list(
         cluster=list(
            title="Cluster",
            title_position="topleft",
            border="black",
            ncol=min(c(length(unique(partition[rows])), 4))
         )
      )
      PHM <- ComplexHeatmap::Heatmap(partition[rows],
         #row_split=partition[rows],
         border=FALSE,
         annotation_legend_param=p_annotation_legend_param,
         use_raster=use_raster,
         col=k_colors,
         name="cluster",
         show_row_names=FALSE,
         width=k_width);
   }

   ## panel_groups
   if (length(panel_groups) > 0) {
      ## Make sure we have some duplicated panel_groups
      if (length(tcount(panel_groups, minCount=2)) > 0) {
         panel_split <- split(seq_along(nmatlist), panel_groups);
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Defining ylims for panel_groups.");
            jamba::printDebug("nmatlist2heatmaps(): ",
               "panel_groups:",
               paste0("(", jamba::cPaste(panel_split), ")"),
               sep="; ")
         }
         ## Calculate shared ylims per panel_group
         panel_ylims <- lapply(panel_split, function(idxs){
            idx_ranges <- lapply(idxs, function(idx){
               nmat <- nmatlist[[idx]][rows,,drop=FALSE];
               itransform <- transform[[idx]];
               if (length(unique(partition)) > 1) {
                  if (verbose>1) {
                     jamba::printDebug("nmatlist2heatmaps(): ",
                        "Calculating colMeans() for panels (",
                        jamba::cPaste(idxs),
                        ") across row clusters.");
                  }
                  plist <- split(names(partition), partition)
                  range(
                     unlist(
                        lapply(plist, function(prows){
                           range(colMeans(
                              itransform(nmat[prows,,drop=FALSE]),
                              na.rm=TRUE))
                        })
                     )
                  )
               } else {
                  if (verbose>1) {
                     jamba::printDebug("nmatlist2heatmaps(): ",
                        "Calculating colMeans() for panels (",
                        jamba::cPaste(idxs),
                        ") across all rows.");
                  }
                  range(
                     colMeans(
                        itransform(nmat),
                        na.rm=TRUE))
               }
            })
            idx_range <- range(pretty(unlist(idx_ranges)));
            idx_range;
         });
         ylims <- panel_ylims[panel_groups];
         if (length(ylims) > 0 && verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "panel_groups ylims: ",
               paste0("(", jamba::cPaste(ylims), ")"),
               sep="; ");
         }
         ## Calculate signal ceiling per panel_group
         panel_ceilings <- lapply(panel_split, function(idxs){
            idx_ceilings <- lapply(idxs, function(idx){
               nmat <- nmatlist[[idx]][rows,,drop=FALSE];
               itransform <- transform[[idx]];
               iceiling <- get_nmat_ceiling(itransform(nmat),
                  signal_ceiling[[idx]],
                  verbose=verbose>1);
            });
            idx_ceiling <- max(unlist(idx_ceilings), na.rm=TRUE);
            if (verbose) {
               jamba::printDebug("nmatlist2heatmaps(): ",
                  "Determined signal ceilings for panels (",
                  jamba::cPaste(idxs),
                  "), ceilings: (",
                  jamba::cPaste(round(digits=2, unlist(idx_ceilings))),
                  "), signal_ceiling=",
                  round(digits=2, idx_ceiling));
            }
            idx_ceiling;
         })
         signal_ceiling <- panel_ceilings[panel_groups];
      }
   }


   ## Iterate each matrix to create heatmaps
   if (length(lens) == 0) {
      lens <- 0;
   }
   lens <- rep(lens, length.out=length(nmatlist));
   if ("gpar" %in% class(axis_name_gp)) {
      axis_name_gp <- rep(list(axis_name_gp),
         length.out=length(nmatlist))
   } else {
      axis_name_gp <- rep(axis_name_gp,
         length.out=length(nmatlist));
   }
   if ("gpar" %in% class(column_title_gp)) {
      column_title_gp <- rep(list(column_title_gp),
         length.out=length(nmatlist));
   } else {
      column_title_gp <- rep(column_title_gp,
         length.out=length(nmatlist));
   }
   if (is.list(axis_name)) {
      axis_name <- rep(axis_name,
         length.out=length(nmatlist));
   } else {
      axis_name <- rep(list(axis_name),
         length.out=length(nmatlist));
   }

   if (length(row_order) == 0) {
      row_order <- TRUE;
   }
   if (is.logical(row_order)) {
      if (any(row_order)) {
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ", sep="",
               c("Defining row_order with ",
                  "EnrichedHeatmap::enriched_score()"));
         }
         row_order <- order(
            EnrichedHeatmap::enriched_score(nmatlist[[main_heatmap]][rows,,drop=FALSE]),
            decreasing=TRUE);
         names(row_order) <- rows;
      } else {
         row_order <- nameVector(rows);
      }
   }
   if (length(row_order) > 1) {
      row_order <- row_order[rows];
   }
   if (any(is.na(row_order))) {
      jamba::printDebug("Fixed NA row_order by assigning rows.");
      row_order <- nameVector(rows);
   }
   #} else if (isFALSE(row_order)) {
   #   row_order <- seq_along(rows);
   #}

   #############################
   ## Iterate each heatmap
   if (verbose) {
      jamba::printDebug("nmatlist2heatmaps(): ",
         "Iterating each heatmap.");
   }
   EH_l <- lapply(seq_along(nmatlist), function(i){
      nmat <- nmatlist[[i]][rows,,drop=FALSE];
      signal_name <- attr(nmat, "signal_name");
      target_name <- attr(nmat, "target_name");
      s_name <- gsub("_at_", "\nat_", signal_name);
      if (length(nmat_names) > 0) {
         signal_name <- nmat_names[i];
      }

      ## For now, do not word wrap, let user do that
      if (1 == 2) {
         signal_name <- gsub("^[\n ]+|[\n ]+$",
            "",
            gsub("\n[\n ]*",
               "\n",
               signal_name))
      }

      color <- nmat_colors[[i]];
      if (length(color) == 0 || is.na(color)) {
         color <- "aquamarine4";
      }
      ## Define ylim
      if (length(ylims) > 0) {
         ylim <- ylims[[i]];
      } else {
         ylim <- NULL;
      }
      itransform <- transform[[i]];
      imat <- itransform(nmat);
      iceiling <- signal_ceiling[[i]];
      divergent <- FALSE;
      if (any(imat < 0)) {
         divergent <- TRUE;
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               c("divergent=", TRUE), sep="")
         }
      }
      if (!is.function(color) && length(color) == 1 && divergent) {
         color2 <- color_complement(color, ...);
         color <- c(color2, middle_color, color);
         divergent <- TRUE;
      }
      if (verbose) {
         if (length(ylim) == 0) {
            ylim_txt <- "NULL";
         } else {
            ylim_txt <- format(ylim, digits=2);
         }
         if (is.function(color)) {
            if (divergent) {
               color_txt <- list(color(c(-10000, 0, 10000)));
            } else {
               color_txt <- list(color(0, 10000));
            }
         } else {
            color_txt <- color;
         }
         jamba::printDebug("nmatlist2heatmaps(): ",
            "signal_name:",
            signal_name,
            ", target_name:",
            target_name,
            ", color:",
            color_txt,
            ", ylim=(", jamba::cPaste(ylim_txt), ")",
            fgText=c(
               rep(list("darkorange",
                  "dodgerblue"),
                  length.out=6),
               NA),
            bgText=c(
               rep(list(NA),
                  length.out=6),
               color_txt)
         );
      }
      if (length(iceiling) > 0 && !is.na(iceiling)) {
         iceiling <- get_nmat_ceiling(imat, iceiling);
         if (divergent) {
            ibreaks <- seq(from=-iceiling,
               to=iceiling,
               length=21);
         } else {
            ibreaks <- seq(from=0,
               to=iceiling,
               length=21);
         }
         if (is.function(color)) {
            colramp <- color;
         } else {
            colramp <- circlize::colorRamp2(
               breaks=ibreaks,
               colors=jamba::getColorRamp(color,
                  defaultBaseColor=middle_color,
                  divergent=divergent,
                  n=21,
                  lens=lens[[i]]));
         }
      } else {
         if (is.function(color)) {
            colramp <- color;
         } else {
            colramp <- jamba::getColorRamp(color,
               n=21,
               defaultBaseColor=middle_color,
               divergent=divergent,
               lens=lens[[i]]);
         }
      }
      if (length(heatmap_legend_param) == 0) {
         legend_width <- grid::unit(3, "cm");
         heatmap_legend_param <- list(direction="horizontal",
            legend_width=legend_width,
            title_position="topleft",
            border="black",
            grid_width=unit(1, "npc"));
      }

      EH <- EnrichedHeatmap::EnrichedHeatmap(imat[rows,],
         split=partition[rows],
         pos_line=FALSE,
         use_raster=use_raster,
         col=colramp,
         border=border[[i]],
         top_annotation=ComplexHeatmap::HeatmapAnnotation(
            lines=EnrichedHeatmap::anno_enriched(
               gp=grid::gpar(col=k_colors),
               value=profile_value,
               ylim=ylim,
               show_error=show_error)
         ),
         heatmap_legend_param=heatmap_legend_param,
         axis_name_gp=axis_name_gp[[i]],
         axis_name=axis_name[[i]],
         axis_name_rot=axis_name_rot,
         name=signal_name,
         column_title=signal_name,
         row_order=row_order[rows],
         column_title_gp=column_title_gp[[i]],
         ...);
      EH;
   });

   ## Optional multi-row layout
   if (hm_nrow > 1 && length(nmatlist) > 1) {
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
               ht_gap=ht_gap,
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
      ht_gap <- rep(ht_gap, length.out=max(c(1, length(nmatlist)-1)));
      if (length(partition) > 0) {
         HM_temp <- PHM + HM_temp;
         main_heatmap_temp <- main_heatmap_temp + 1;
         ht_gap <- grid::unit.c(grid::unit(1, "mm"), ht_gap);
      }
      if (length(AHM) > 0) {
         HM_temp <- AHM + HM_temp;
         main_heatmap_temp <- main_heatmap_temp + 1;
         ht_gap <- grid::unit.c(grid::unit(1, "mm"), ht_gap);
      }
      if (length(MHM) > 0) {
         HM_temp <- HM_temp + MHM;
         ht_gap <- grid::unit.c(ht_gap, grid::unit(1, "mm"));
      }
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "ht_gap:");
         print(ht_gap);
      }
      if (length(title) > 0 || length(caption) > 0) {
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
      } else {
         ComplexHeatmap::draw(HM_temp,
            ht_gap=ht_gap,
            main_heatmap=main_heatmap_temp);
      }
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
#' @family jam import functions
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
#' Recognized transform names:
#'
#' * `"none"` or `"linear"` returns the data without change
#' * `"log2signed"` applies `jamba::log2signed()`, defined as
#'      `log2(1+x)` transform to the absolute value, then multiplies
#'      by the original `sign(x)`
#' * `"exp2signed"` applies the inverse of `"log2signed"`, which
#'    exponentiates numeric values which were previously transformed
#'    with `log2(1+x)`. Note that value=0 when exponentiated becomes `+1`
#' * `"sqrt"` applies square root transform, equivalent to `sqrt()`
#' * `"cubert"` applies cube root `x^(1/3)`
#' * `"qrt"` applies fourth root `x^(1/4)`
#' * `"frt"` or `"fthrt"` applies fifth root `x^(1/5)`
#' * `"square"` applies `x^2` to absolute value, multiplied by the `sign(x)`
#' * `"cube"` applies `x^3`;
#'
#' Any other character name is used to find a function with the same
#' name, and if the function is found it is used. The function is
#' expected to take input `x` and return corresponding output with the
#' same length as the input. A function such as `max()` does not fit
#' this criteria, but a function such as `log2()` is acceptable.
#'
#' @return `function` or `NULL` when no matching function is
#'    found, or `list` is returned when the input `transform`
#'    has multiple values.
#'
#' @family jam utility functions
#'
#' @param transform `character` string or `function`, or `list`
#'    that may contain `character` string or `function`.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' get_numeric_transform("log2signed")
#'
#' transform_list <- get_numeric_transform(
#'    c("none",
#'    "log2signed",
#'    log2,
#'    sqrt));
#' jamba::sdim(transform_list);
#'
#' x <- jamba::nameVector(0:16);
#' transform_list[[1]](x)
#'
#' transform_list[[2]](x)
#'
#' @export
get_numeric_transform <- function
(transform,
 ...)
{
   #
   if (length(transform) > 1) {
      if (length(names(transform)) == 0) {
         ## Determine names from the input functions
         sys <- sys.call();
         sys2 <- as.list(as.list(sys)[[2]])[-1];
         sys3 <- lapply(sys2, function(i){
            if ("call" %in% class(i)) {
               head(deparse(i), 1)
            } else {
               as.character(i)
            }
         })
         trans_names <- jamba::cPaste(sys3)
         names(transform) <- trans_names;
      }
      transform <- lapply(transform, function(it){
         get_numeric_transform(it, ...)
      });
      return(transform);
   }
   it <- transform;
   if (is.atomic(it)) {
      if ("none" %in% it) {
         it <- function(x)x;
      } else if ("log2signed" %in% it) {
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
      } else if (any(c("qrt", "quadrt") %in% it)) {
         it <- function(x){sign(x)*abs(x)^(1/4)};
      } else if (any(c("frt", "fthrt") %in% it)) {
         it <- function(x){sign(x)*abs(x)^(1/5)};
      } else {
         itf <- tryCatch({
            get(it);
         }, error=function(e){
            jamba::printDebug("get_numeric_transform(): ",
               "Error:",
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


#' Helper function to calculate signal ceiling of numeric matrix
#'
#' Helper function to calculate signal ceiling of numeric matrix
#'
#' This function is called by `nmatlist2heatmaps()` and is not
#' intended to be called directly.
#'
#' @family jam utility functions
#'
#' @export
get_nmat_ceiling <- function
(imat,
 iceiling=NULL,
 verbose=TRUE,
 ...)
{
   if (length(iceiling) == 0 || any(is.na(iceiling))) {
      iceiling <- max(abs(imat),
         na.rm=TRUE);
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "   Applied max(nmat) ceiling=",
            round(digits=3, iceiling));
      }
   } else if (iceiling > 0 && iceiling <= 1) {
      # apply quantile
      imat_values <- setdiff(abs(imat), 0);
      iquantile <- quantile(imat_values,
         probs=iceiling,
         na.rm=TRUE);
      if (verbose) {
         jamba::printDebug("get_nmat_ceiling(): ",
            "Applied quantile=",
            iceiling,
            " and defined ceiling=",
            round(digits=3, iquantile));
      }
      iceiling <- iquantile;
   } else if (verbose) {
      jamba::printDebug("nmatlist2heatmaps(): ",
         "   Applied ceiling=",
         round(digits=3, iceiling));
   }
   return(iceiling);
}

#' Create color complement by rotating the color hue
#'
#' Create color complement by rotating the color hue
#'
#' This function rotates the color hue to create a complementary
#' color for each `color` input. It differs from standard methods
#' by using warped color hue by default (`useWarpHue=TRUE`), which
#' uses a red-yellow-blue color wheel instead of R default
#' red-green-blue. It also imposes a minimum chroma, which
#' ensures the output color is reasonably high in color
#' saturation.
#'
#' @family jam utility functions
#'
#' @param color `character` vector of R compatible colors.
#' @param Hflip numeric value in degrees (from 0 to 360) added
#'    to the color hue to produce the final color hue.
#' @param Cfloor numeric value used to limit output chroma `C`
#'    values to this minimum value.
#' @param Lrange `numeric` vector with the allowed range of output
#'    luminance `L` values. When supplied, output values are
#'    simply forced to this range with no other scaling of intermediate
#'    values.
#' @param useWarpHue `logical` indicating whether to use the warp
#'    hue functions `colorjam::h2hw()` and `colorjam::hw2h()` which
#'    effectively change the color wheel from red-green-blue to
#'    red-yellow-blue.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' n <- 5;
#' rc <- colorjam::rainbowJam(n);
#' rc_comp <- color_complement(rc, Cfloor=180);
#' rc_comp2 <- color_complement(rc, Cfloor=180, useWarpHue=FALSE);
#' jamba::showColors(list(rainbowJam=rc,
#'    `complement_warped\n(default uses RYB)`=rc_comp,
#'    `complement_nowarp\n(unwarped uses RGB)`=rc_comp2));
#'
#' ## divergent color gradients through white
#' ## hint: use higher lens value to make middle colors more intense
#' rc_ramps <- lapply(jamba::nameVector(seq_along(rc)), function(i){
#'    j <- jamba::getColorRamp(c(rc[i], "white", rc_comp[i]),
#'       n=25,
#'       lens=0,
#'       divergent=TRUE);
#'    names(j) <- "";
#'    names(j)[1] <- "original colors";
#'    names(j)[25] <- "color complements";
#'    j;
#' });
#' jamba::showColors(rc_ramps, groupCellnotes=TRUE, groupByColors=FALSE);
#'
#' ## divergent color gradients through white
#' ## hint: use higher lens value to make middle colors more intense
#' rc_ramps2 <- lapply(jamba::nameVector(seq_along(rc)), function(i){
#'    j <- jamba::getColorRamp(c(rc[i], "black", rc_comp[i]),
#'       n=25,
#'       lens=1,
#'       divergent=TRUE);
#'    names(j) <- "";
#'    names(j)[1] <- "original colors";
#'    names(j)[25] <- "color complements";
#'    j;
#' });
#' jamba::showColors(rc_ramps2, groupCellnotes=TRUE, groupByColors=FALSE);
#'
#' @export
color_complement <- function
(color,
 Hflip=180,
 Cfloor=160,
 Lrange=c(0, 100),
 useWarpHue=TRUE,
 ...)
{
   hcl <- jamba::col2hcl(color);
   H <- hcl["H",];
   if (useWarpHue) {
      H <- colorjam::h2hw(H);
   }
   newH <- (Hflip + H) %% 360;
   if (useWarpHue) {
      newH <- colorjam::hw2h(newH);
   }
   hcl["H",] <- newH;
   if (length(Cfloor) > 0) {
      hcl["C",] <- jamba::noiseFloor(hcl["C",],
         minimum=Cfloor);
   }
   if (length(Lrange) > 0) {
      hcl["L",] <- jamba::noiseFloor(hcl["L",],
         minimum=min(Lrange),
         ceiling=max(Lrange));
   }
   color2 <- jamba::hcl2col(hcl);
   return(color2);
}
