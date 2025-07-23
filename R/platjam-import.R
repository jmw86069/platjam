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
#' @family jam coverage heatmap functions
#' @family jam import functions
#'
#' @returns `normalizedMatrix` numeric matrix, where additiona
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
#'       top_annotation=ComplexHeatmap::HeatmapAnnotation(
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
#'       width=grid::unit(3, "mm"));
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
         xseq <- jamba::nameVector(seq_along(x), names(x));
         nmatlist <- lapply(xseq, function(i){
            if (length(signal_name) >= i) {
               use_name <- signal_name[i];
            } else {
               use_name <- names(x)[i]
            }
            coverage_matrix2nmat(x=x[i],
               filename=NULL,
               signal_name=use_name,
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
         stop("Must supply x as a data.frame, or supply one filename.");
      } else if (!file.exists(head(x, 1))) {
         stop("No x supplied as data.frame, and filename is not found.");
      } else {
         if (verbose) {
            jamba::printDebug("coverage_matrix2nmat(): ",
               "Importing data from filename:",
               filename);
         }
         x <- data.table::fread(head(filename, 1),
            sep="\t",
            data.table=FALSE);
         if (length(signal_name) == 0) {
            if (length(names(filename)) > 0 && nchar(head(names(filename), 1)) > 0) {
               signal_name <- head(names(filename), 1);
            } else {
               signal_name <- gsub("[.](matrix|matrix.gz)$", "",
                  basename(head(filename, 1)));
            }
         }
      }
   } else if (jamba::igrepHas("character", class(x)) &&
         file.exists(head(x, 1))) {
      if (verbose) {
         jamba::printDebug("coverage_matrix2nmat(): ",
            "Importing data from filename sent as x:",
            x);
      }
      filename <- x;
      x <- data.table::fread(head(filename, 1),
         sep="\t",
         data.table=FALSE);
      if (length(signal_name) == 0) {
         if (length(names(filename)) > 0 && nchar(head(names(filename), 1)) > 0) {
            signal_name <- head(names(filename), 1);
         } else {
            signal_name <- gsub("[.](matrix|matrix.gz)$", "",
               basename(head(filename, 1)));
         }
      }
   }
   if (verbose) {
      jamba::printDebug("coverage_matrix2nmat(): ",
         "signal_name:",
         signal_name);
      jamba::printDebug("coverage_matrix2nmat(): ",
         "names(filename):",
         names(filename));
   }

   if (length(signal_name) == 0) {
      signal_name <- "signal";
   }

   if (!jamba::igrepHas("data.frame|dataframe|tibble|data.table", class(x))) {
      stop(paste0("Supplied x must be class ",
         "data.frame, DataFrame, data.table, or tibble."));
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
      jamba::printDebug("coverage_matrix2nmat(): ",
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
#' # Annotation Data
#'
#' When `anno_df` is provided as a `data.frame` the rows are synchronized
#' alongside the heatmap rows. Column values are color-coded, categorical
#' for `character` columns, and using color gradient for `numeric` columns.
#'
#' Rows can optionally be split by argument `partition`, which can be a vector
#' of group values associated with rows, or one or more columns in
#' `colnames(anno_df)` whose values are used to sub-divide the rows.
#'
#' # Row Clustering / Partitioning
#'
#' Rows can be clustered using k-means clustering with argument `k_clusters`.
#' By default it uses `k_method="correlation"`, which applies a novel
#' and effective correlation metric, clustering row data by the profile shape.
#' The typical default, which is used when the `amap` R package is not
#' installed, is to use `"euclidean"` distance, which tends to cluster
#' based upon signal magnitude moreso than the shape.
#'
#' When k-means clustering `k_clusters` and `partition` are both enabled,
#' each partition is independently k-means clustered, which improves
#' results compared to applying global k-means before applying partitions.
#' Use `min_rows_per_k` to adjust the relative number of `k` clusters
#' based upon the number of observed rows.
#'
#' # Display Layout
#'
#' Heatmaps are arranged in the following order, dependent upon
#' the data provided:
#'
#' * Annotation heatmap, if `anno_df` is provided.
#'
#'    * Color assignment can be provided using `color_sub` either as
#'    a named vector of R colors whose names match values in each column,
#'    or as a `list` named by `colnames(anno_df)`, with named color
#'    assignments, or a color `function` for `numeric` columns.
#'
#' * Partition heatmap, if `partition` is provided.
#' * Enrichment heatmaps, one for each entry in `nmatlist`.
#'
#'    * Above each heatmap is the metaplot, drawn using
#'    `EnrichedHeatmap::anno_enriched()`.
#'    * When `partition` and/or `k_clusters` are defined,
#'    the plot will include one profile line for each row grouping.
#'    * When `show_error=TRUE` each line will also be shaded using 95%
#'    standard deviation.
#'    * The heatmap color gradient is applied starting at zero, extending to
#'    `signal_ceiling` for each heatmap. When `signal_ceiling` is <=1 it
#'    uses the quantile of non-zero values in the matrix data, otherwise
#'    it applies a fixed numeric maximum. Numeric values above the
#'    `signal_ceiling` threshold are colored using the maximum color.
#'    * When there are negative values, the color key uses a divergent
#'    color scale. When `nmat_colors` value for the heatmap is a single color,
#'    the complementary color is used for negative values; otherwise it is
#'    assumed to define a divergent color scale.
#'    * The y-axis range on metaplots is defined by observed values, and
#'    when `panel_groups` is defined, the y-axis `ylim` is shared among
#'    all heatmaps in each panel group.
#'
#' * Marked row heatmap, if `anno_row_marks` is provided. It uses an empty
#' heatmap, associated with row mark annotations for a subset of row labels,
#' in the same order as the coverage heatmaps.
#' * Color legends are displayed in the same order:
#'
#'    * annotation colors for each column in `anno_df`
#'    * partition/cluster colors
#'    * color gradients for each coverage heatmap in order, or
#'    when `panel_groups` is provided it displays the color key for
#'    only the first heatmap in each panel group.
#'
#' @returns `list` with heatmap components that can be reviewed, or
#'    optionally rendered into a figure:
#'    * `"AHM"`: annotation heatmap, when `anno_df` is supplied
#'    * `"PHM"`: partition heatmap, when partitioning and/or k-means clustering
#'    is used
#'    * `"EH_l"`: `list` of `ComplexHeatmap::Heatmap` objects
#'    * `"MHM"`: marked heatmap, containing optional row labels
#'    * `"HM_drawn"`: when `hm_nrow=1` this is the output after drawing the
#'    heatmap, in the form: `ComplexHeatmap::HeatmapList`. This object can
#'    be drawn again if needed, or used to determine exact row orders.
#'    * `"fn_params"`: `list` of useful function parameters, including
#'    some calculated during processing such as `panel_groups`, `ylims`,
#'    `signal_ceiling`, etc.
#'    * `"hm_caption"`: `character` version of heatmap captions
#'    * `"adjust_df"`: `data.frame` when `recenter_heatmap` or
#'    `restrand_heatmap` are defined, which contains a summary of each
#'    row, with colnames:
#'    `"summit_name"` for recentering; and
#'    `"restrand"` for restranding.
#'
#' @param nmatlist `list` containing `normalizedMatrix` objects,
#'    usually the output from `coverage_matrix2nmat()`.
#' @param panel_groups `character` vector with values for each `nmatlist`
#'    entry, which defines groups of heatmap panels.
#'    Each panel group shares:
#'    * numeric range for the heatmap color gradient, defined by the first
#'    `signal_ceiling` value for the group. Standard rules apply, such that
#'    values below 1 represent a quantile signal threshold, and values above 1
#'    represent a fixed numeric threshold.
#'    * one color key, labeled by `names(panel_groups)` to represent all
#'    panels in the group
#'    * `ylim` y-axis range for the profile plot, either determined dynamically
#'    or by the first `ylim` provided for the panel group
#'    * When `nmat_colors` is not defined, each panel group is assigned
#'    one categorical color which is applied to all heatmaps in the group.
#'    * When `nmat_colors` is defined, each panel uses the color as defined,
#'    however the color key only uses the color gradient from the first
#'    panel in the group.
#' @param title,caption `character` string used as an overall title or
#'    caption, respectively, displayed at the top of all heatmap output.
#' @param title_gp `grid::gpar` object to customize the title fontsize,
#'    fontface, color (col), etc.
#' @param upstream_length,downstream_length `numeric` (optional) range of
#'    coordinates to display across all heatmaps. This argument is intended
#'    when the input `nmatlist` contains a wider range of coordinates
#'    than should be displayed. The columns in `nmatlist` are subset
#'    to retain only those columns within the range `downstream_length`
#'    to `upstream_length`, assuming the middle coordinate is zero.
#'    This step calls `zoom_nmatlist()`.
#'    Note this step does not expand the displayed region.
#' @param k_clusters `integer` number of k-means clusters to
#'    use to partition each heatmap. Use `0` or `NULL` for
#'    no clustering (default).
#'    Note `k_clusters` can be a `numeric` vector, in which case it is applied
#'    across unique groups defined by `partition` if provided.
#'    If `names(k_clusters)` match values in `partition` they will be applied
#'    by name, otherwise they are applied in the order the clusters are
#'    defined by `partition`.
#'    Each group is clustered to that many k clusters, provided it
#'    also meets the threshold `min_row_per_k` - which is intended to prevent
#'    clustering 10 rows into 10 k-means clusters.
#' @param min_rows_per_k `numeric` minimum rows required per k-means
#'    cluster, used only when `k_clusters` is greater than 1.
#'    With default `min_rows_per_k=10`, a partition with 100 or fewer rows
#'    can only have `k=1`, and partition with 101 rows can have `k=2`.
#'    This limit protects from k-means clustering small partitions.
#' @param k_subset `integer` vector of k-means clusters to retain.
#'    This argument is intended to "zoom in" (or "drill down") to one
#'    or more k-means clusters of interest.
#'    When both `k_clusters` and `partition` are provided, this argument
#'    must exactly match the row title as displayed in the heatmap.
#' @param k_colors `character` vector of R colors, or `NULL` to use
#'    the output of `colorjam::rainbowJam(k_clusters)`.
#'    These colors are applied to `k_clusters` and/or `partition`:
#'    * When `partition` is provided, `names(k_colors)` are used when present,
#'    otherwise colors are assigned in order of `partition` groups.
#'    When `k_clusters` is also defined, each partition color is split into
#'    a light-to-dark gradient based upon the number of k_clusters.
#'    * When `partition` is not provided, `k_colors` are applied to k-means
#'    clusters in the order the colors are provided.
#' @param k_width `unit` width of the k-means cluster color
#'    bar, used with `k_clusters`, default is 5 mm width.
#' @param k_method `character` string indicating the distance
#'    used by k-means, where the common default is
#'    `"euclidean"`, however a useful alternative for
#'    sequence coverage data is `"correlation"` as implemented
#'    in `amap::Kmeans()`. Available methods:
#'    * `"euclidean"` (default) calculates the typical Euclidean distance,
#'    which tends to emphasize total signal moreso than the specific
#'    shape of the signal.
#'    * `"correlation"` when the R package `amap` is available, this method
#'    emphasizes the shape of signal profiles, and is particularly effective.
#'    It is also called "centered Pearson" since data is centered prior
#'    to calculating correlation.
#'    * `"pearson"` when the R package `amap` is available, this method
#'    is also called "not centered Pearson" since data is not centered
#'    prior to calculating correlation.
#'    * `"spearman"` when the R package `amap` is available, this method
#'    computes distance based upon rank differences. It has not been tested
#'    much in this context.
#' @param k_heatmap `integer` with one or more values indicating which
#'    `nmatlist` entries to use for k-means clustering,
#'    default uses `main_heatmap`. This value is only used when
#'    `k_clusters` is greater than 1. This argument is useful for
#'    clustering multiple coverage heatmaps together.
#' @param partition `character` or `factor` vector used to split rows
#'    of each matrix in `nmatlist`, and **must named by rownames**
#'    in `nmatlist`.
#'    This value is converted to `factor`, and will honor provided
#'    factor levels if already defined.
#'    * When `partition` and `k_clusters` are both defined, the
#'    data is first grouped by `partition` then each partition group
#'    is separately k-means clustered, using rules described for
#'    `k_clusters` and `min_rows_per_k`.
#'    Colors from `k_colors` are assigned to each partition value,
#'    then colors are split to light-to-dark gradient
#'    using `jamba::color2gradient()`.
#' @param row_title_rot `numeric` value in degrees, to rotate the
#'    partition labels on the left, when either `partition` or
#'    `k_clusters` are provided. The default `0` uses horizontal
#'    text. For long labels, it may be better to use `30` or `60`.
#' @param partition_counts `logical` indicating whether to include the
#'    number of rows in each partition, default `TRUE`.
#'    Note that this setting is active if `k_clusters` and/or `partition`
#'    are supplied. Any situation where rows are split, the number of
#'    rows will be displayed.
#' @param partition_count_template `character` format used when
#'    `partition_counts=TRUE`, used together with `glue::glue()` to format
#'    each row partition. The default: `"{partition_name}\n({counts} rows)"`
#'    will print for example: `"A\n(125 rows)"`
#' @param rows optional vector to define subset rows, or specific row order:
#'    * `character` vector of rownames in `nmatlist`, or
#'    * `integer` vector with row numbers (row index) values.
#'
#'    Note that even when using a subset of `rows` the data may also be
#'    subset based upon available `names(partition)` and
#'    `rownames(anno_df)`.
#' @param row_order `integer` vector used to order rows, intended to
#'    allow ordering data based upon a specific heatmap, or using
#'    different logic than the default.
#'    * When `row_order=NULL` (default) or `row_order=TRUE` it calls
#'    `EnrichedHeatmap::enriched_score()` using data from `main_heatmap`.
#'    When there are multiple values for `main_heatmap` (which is default),
#'    then scores are calculated for each matrix, then the average score
#'    is used per row.
#'
#'    The `enriched_score()` function generates a weighted score with
#'    heighest weight at the center position, with progressively lower
#'    weight working outward where the maximum distance has zero weight.
#'    The technique sorts signal which emphasizes highest enriched signal
#'    at the center of the matrix.
#'    * When `row_order=FALSE` the data is ordered in the same order
#'    they appear in `nmatlist`, or when `anno_df` and `byCols` are supplied,
#'    the rows in `anno_df` are sorted using
#'    `jamba::mixedSort(anno_df, byCols=byCols)` and the resulting row order
#'    is used.
#' @param nmat_colors `character` vector of R colors,
#'    to colorize each heatmap.
#'    * When `nmat_colors=NULL` (default) and `panel_groups` is not defined,
#'    `colorjam::rainbowJam()` is used to assign one unique color
#'    to each heatmap panel.
#'    * When `nmat_colors=NULL` and `panel_groups` is defined,
#'    `colorjam::rainbowJam()` is used to assign one unique color
#'    to each unique panel group, and the same color is applied to each
#'    heatmap panel in each panel group.
#' @param middle_color `character` R color, default `middle_color="white"`,
#'    used as the middle color when creating a divergent color gradient.
#'    This color should usually be either `"white"` or `"black"`, but
#'    sometimes can be slightly off-white or off-black to apply some
#'    distinction from the background color.
#' @param nmat_names `character` vector, or `NULL`, optional,
#'    used as custom names for each heatmap in `nmatlist`.
#'    When `nmat_names=NULL` the `signal_name` values are
#'    used from each `nmatlist` entry attribute: `attr(nmat, "signal_name")`
#' @param main_heatmap `integer` index to define one or more entries
#'    in `nmatlist` as the main heatmap used for clustering and row ordering.
#'    Note that `k_heatmap` will override this option when provided.
#'    By default `main_heatmap=NULL` will cause all heatmaps to be used for
#'    row ordering.
#' @param anno_df `data.frame` or object that can be coerced to `data.frame`
#'    whose `rownames(anno_df)` must match rownames in the nmatlist data.
#'    When `rownames(anno_df)` does not match, this function fails with
#'    an error message.
#'    * Data can optionally be sorted by defining `byCols`.
#'    * When provided, data in `nmatlist` is automatically subsetted
#'    to the matching `rownames(anno_df)` also present in `nmatlist`.
#'    * When `rows` is also defined, the data will be subsetted by
#'    the `rows` and by the `rownames(anno_df)` present in `nmatlist`.
#' @param byCols `character` vector of `colnames(anno_df)` used to
#'    sort the `data.frame`. This argument is passed to
#'    `jamba::mixedSortDF()` and follows its rules, for example prefix `"-"`
#'    causes the column to be sorted in reverse. Multiple columns can be
#'    sorted, in the order they are provided, and factor levels are
#'    honored for factor columns.
#' @param color_sub accepts input in two forms:
#'    1. `character` vector of R colors named by `character` values
#'    2. `list` output from `design2colors()` where each `list` element
#'    is named by colnames present in `anno_df`,
#'    and each `list` value is either:
#'
#'       * `character` vector of colors named by `character` value, or
#'       * color `function` as defined by `circlize::colorRamp2()`,
#'       which takes a `numeric` value and returns a `character` R color.
#'
#'    * When values for any column in `anno_df` does not have colors
#'    assigned by one mechanism above, colors are assigned using
#'    `colorjam::group2colors()`.
#'    * When `partition` is defined, colors are assigned either by
#'    matching unique partition values with `names(color_sub)`,
#'    or with `attr(color_sub, "color_sub")` if present, which may contain
#'    the full set of name-color assignments when `color_sub` is
#'    provided as a `list`. Otherwise if `color_sub` is provided as a `list`
#'    each entry is compared with `partition` values until values can
#'    be fully matched. Failing these steps, colors are assigned to
#'    unique `partition` values, then if `k_clusters` is also supplied,
#'    the partition colors are then split by `colorjam::color2gradient()`
#'    across the k-means clusters for each partition.
#' @param anno_row_marks `character` optional vector of `rownames`
#'    in `nmatlist` that should be labeled beside the heatmaps using
#'    `ComplexHeatmap::anno_mark()`.
#'    * Note `anno_row_labels` can be used to supply custom labels,
#'    or one or more columns in `anno_df`.
#'    * When `anno_row_labels=NULL` (default) it displays the value
#'    in `anno_row_marks` itself.
#' @param anno_row_labels `character` vector of optional labels to use
#'    when `anno_row_marks` is supplied.
#'    * When `anno_row_labels=NULL` (default) it uses rownames defined
#'    in `anno_row_marks.
#'    * It can be a `character` vector of actual labels, with names
#'    that match `anno_row_marks` (thus rownames in `nmatlist`).
#'    * It can be a `character` vector with one or more `colnames(anno_df)`,
#'    which creates labels by concatenating values across columns,
#'    delimited with space `" "`.
#' @param anno_row_gp `grid::gpar` object used to customize the text label
#'    displayed when `anno_row_marks` is defined. The default fontsize 14
#'    is intended to be larger than other default values, for legibility.
#' @param recenter_heatmap,recenter_range,recenter_invert arguments
#'    are passed to `recenter_nmatlist()` to apply re-centering.
#'    * Note that recenter will always occur before restrand.
#' @param summit_names `character` default NULL, optional colnames to
#'    use for recentering, which applies a previously defined set of
#'    summit positions to use. It ignores all other recenter arguments.
#' @param restrand_heatmap,restrand_range,restrand_buffer,restrand_invert
#'    arguments are passed to `restrand_nmatlist()` to apply re-stranding.
#'    * Note that recenter will always occur before restrand.
#' @param top_annotation `HeatmapAnnotation` or `logical` or `list`:
#'    * `top_annotation=TRUE` (default) uses the default
#'    `EnrichedHeatmap::anno_enriched()` to display the signal profile
#'    for each row partition and/or k-means cluster.
#'    * `top_annotation=FALSE` does not display a top annotation.
#'    * object `HeatmapAnnotation` as produced by
#'    `ComplexHeatmap::HeatmapAnnotation(EnrichedHeatmap::anno_enriched())`
#'    or equivalent. This form is required for the annotation
#'    function to be called successfully on each heatmap in `nmatlist`.
#'    * a `list` of objects to be applied sequentially to
#'    each `nmatlist` coverage heatmap in order, intended to allow custom
#'    top annotation for each heatmap.
#' @param top_anno_height `unit` object to define the default
#'    height of the `top_annotation`. When `top_annotation`
#'    is not defined, the default method uses
#'    `EnrichedHeatmap::anno_enriched()` with
#'    `height=top_anno_height`.
#' @param top_axis_side `character` value indicating which side
#'    of the top annotation to place the y-axis labels.
#'    * When only one value is defined, it is recycled across `nmatlist`.
#'    * Otherwise it is used when `panel_groups` are defined,
#'    and the top annotation is labeled for only one panel in
#'    each panel group using the side as defined. Labels are displayed
#'    for each contiguous set of panel groups, so that heatmaps
#'    in the same panel group can be ordered in different subsets.
#'    Consider panel groups in this order: A, A, B, B, A, A. It would
#'    display one set of axis labels for the first two panels in A, then
#'    one axis label for the next two panels in B, then one axis label
#'    again for the final two panels in A.
#'    * Values should be one of:
#'
#'       * `"left"`,`"right"`: axis labels on this side of each panel group
#'       * `"both"`: axis labels on both sides of each panel group, useful
#'       when panel groups have a fairly large number of panels.
#'       * `"none"`: display no axis labels
#'       * `"all"`: display axis labels for every panel even within panel group.
#' @param legend_max_ncol `integer` number indicating the maximum
#'    number of columns allowed for a categorical color legend.
#' @param legend_base_nrow `integer` number indicating the base
#'    number of rows used for a categorical color legend, before
#'    additional columns are added. Once the number of elements
#'    exceeds `(legend_max_ncol * legend_base_nrow)` then
#'    rows are added, but columns never exceed `legend_max_ncol`.
#' @param legend_max_labels `integer` to define the maximum labels
#'    to display as a color legend. When any `anno_df` column contains
#'    more than this number of categorical colors, the legend is
#'    not displayed, in order to prevent the color legend from filling
#'    the entire plot device, thus hiding the heatmaps.
#' @param show_heatmap_legend `logical` indicating whether to display the
#'    color legend for each heatmap entry in `nmatlist`. When `panel_groups`
#'    are supplied, color legends are displayed only for the first
#'    heatmap in each unique panel group, unless `show_heatmap_legend=FALSE`,
#'    or unless `show_heatmap_legend` is already defined for every heatmap.
#' @param heatmap_legend_param `list` with optional heatmap legend settings.
#'    By default `NULL` causes this argument to be defined internally,
#'    however when provided it overrides any internal settings and is used
#'    directly. The `list` should be `length(nmatlist)`, or is recycled
#'    to that length.
#' @param heatmap_legend_direction `character` string used when
#'    `show_heatmap_legend=TRUE` and `heatmap_legend_param` is not already
#'    provided.
#'    * By default `heatmap_legend_direction="horizontal"` displays
#'    the color gradient in the legend horizontally as a continuous scale,
#'    with labels defined in `EnrichedHeatmap::EnrichedHeatmap()`, and
#'    width equal to `grid::unit(1, "npc")` which uses the full width of
#'    the color legend area.
#'    * When `heatmap_legend_direction="vertical"` the color legend is
#'    displayed vertically, with width `grid::unit(5, "mm")`.
#' @param annotation_legend_param `list` optional parameters passed to
#'    the annotation legend functions, intended to provide customization.
#'    The `list` should be named by each annotation entry to be customized,
#'    and any annotation entries not defined `annotation_legend_param`
#'    use the default behavior of `ComplexHeatmap::HeatmapAnnotation()`,
#'    which will assign its own set of colors and use default legend
#'    parameters by default.
#'    When `annotation_legend_param=NULL` (default) then all colors
#'    are defined, and all legends are displayed using this function
#'    defaults. When there are more labels than `legend_max_labels`
#'    the color legend will be hidden for that annotation legend entry.
#' @param hm_nrow `integer` number of rows used to display
#'    the heatmap panels. This mechanism is somewhat experimental,
#'    and is used to split a large number of coverage heatmaps into
#'    two rows of heatmaps.
#'    * The matrix data row order is consistent across all heatmap panels.
#'    * The annotation data is displayed to the left of each row of
#'    heatmap panels.
#' @param transform one of the following:
#'    * `character` string referring to a numeric transformation,
#'    passed to `get_numeric_transform()`. Commonly used strings:
#'
#'       * `"log2signed"` calls `jamba::log2signed()`, which applies
#'       `log2(1+x)` to the absolute value, multiplied by `sign(x)`
#'       * `"sqrt"` applies square root to the absolute value, multiplied
#'       by the `sign(x)`
#'       * `"cubert"` applies cube root `x^(1/3)`
#'       * `"qrt"` applies fourth root `x^(1/4)` to the absolute value,
#'       multiplied by the `sign(x)`
#'
#'    * `function` that applies a numeric transformation.
#'    Valid `character` string values:
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
#' @param transform_label `character` optional vector of transformation labels
#'    to use. When `transform_label=NULL` (default) it uses `names(transform)`
#'    if present, then the `character` string of `transform`, otherwise
#'    is left blank. When `transform="none"` no label is displayed.
#'    By default, transform labels are surrounded by parentheses, for example
#'    `"(log2signed)"` and placed on a new line below each coverage heatmap
#'    title. To suppress the transformation in the title, supply
#'    `transform_label=""`.
#' @param signal_ceiling `numeric` vector whose values are recycled
#'    to length `length(nmatlist)`. The `signal_ceiling`
#'    defines the maximum numeric value to the color ramp for
#'    each matrix in `nmatlist`. The value is passed to `get_nmat_ceiling()`,
#'    which recognizes three numeric forms:
#'    1. `signal_ceiling=NULL`: (default) the maximum absolute value
#'    is used as the ceiling.
#'    2. `signal_ceiling > 1`: the specific numeric value
#'    is applied as a fixed ceiling, even if the value is above or below
#'    the maximum absolute value in the data matrix. This setting is useful
#'    for defining a fixed meaningful threshold across `nmatlist` entries.
#'    3. `signal_ceiling > 0` and `signal_ceiling <= 1`: the numeric value
#'    defines a quantile threshold calculated using signal in the data matrix,
#'    excluding values of zero. For example `signal_ceiling=0.75` calculates
#'    ceiling `quantile(x, probs=0.75)`, using non-zero values.
#'
#'    Note that the ceiling is only applied to the color scale and
#'    not to the underlying data. The row clustering and row ordering
#'    steps use the full data range, after applying the appropriate
#'    `transform` where applicable.
#'
#'    To apply a numeric ceiling to the data itself, it should be done
#'    at the level of `nmatlist` beforehand.
#' @param lens `numeric` adjustment to the intensity of the color gradient,
#'    used only when the corresponding `nmat_colors` entry uses a fixed
#'    set of colors. `lens` above zero create more rapid color changes,
#'    making the gradient more visually intense, values below zero reduce
#'    the intensity.
#'    The `lens` values are recycled to `length(nmatlist)` as needed.
#'    Note that `signal_ceiling` defines the `numeric` value at which
#'    the maximum color is applied, while `lens` adjusts the intensity of
#'    the intermediate values in the color gradient.
#' @param anno_lens `numeric` value used to scale the annotation
#'    heatmap color scales, see `lens` for details. This value is applied
#'    to `numeric` columns only when `anno_df` is provided.
#' @param axis_name `character` string with optional custom label used for
#'    the target region label in each heatmap panel.
#'    * When `axis_name=NULL` (default), the `attr(nmat, "target_name")`
#'    label will be used, which is usually "target", along with the
#'    upstream and downstream length as stored in `attr(nmat, "extend")`.
#'    * a `character` vector will be applied as the center
#'    (target) label on each heatmap, using the upstream and downstream
#'    length as stored in `attr(nmat, "extend")`.
#'    * a `list` is expected to have three labels per vector element,
#'    corresponding to the upstream, target, and downstream axis label.
#'    This `list` is recycled to `length(nmatlist)`.
#' @param axis_name_gp object of `grid::gpar` applied to the x-axis label
#'    graphic parameters. For example, to customize the x-axis font size,
#'    use the form: `grid::gpar(fontsize=8)`.
#' @param axis_name_rot `numeric` value either `0` or `90` indicating
#'    whether to rotate the x-axis names below each heatmap, where
#'    `axis_name_row=90` (default) will rotate labels vertically,
#'    and `axis_name_row=0` will display labels horizontally.
#'    * Note that `axis_name_rot` also controls the rotation of
#'    annotation (`anno_df`) and partition (`partition` or `k_clusters`)
#'    annotation labels, below each annotation heatmap.
#' @param column_title_gp object `grid::gpar` or `list` of `grid::gpar`
#'    objects, applied across entries in `nmatlist` to customize the title
#'    displayed above each heatmap panel.
#'    For example to alter the font size, use `grid::gpar(fontsize=14)`.
#'    This argument is passed to `ComplexHeatmap::Heatmap()`, and can
#'    be customized for each heatmap as needed.
#' @param seed `numeric` value used with `set.seed()` to
#'    set the random seed. Set to `NULL` to avoid running
#'    `set.seed()`.
#' @param ht_gap `unit` size to specify the gap between multiple heatmaps.
#'    This argument is passed to `ComplexHeatmap::draw()`. An example
#'    is `grid::unit(8, "mm")` to specify 8 millimeters.
#' @param row_anno_padding,column_anno_padding,legend_padding `grid::unit`
#'    to define the padding between heatmap body, and row annotation,
#'    column annotation, and heatmap color legend, respectively.
#'    * The default values are intended to provide more space between heatmap
#'    and these features than between heatmap subsections (`row_gap`).
#'    The defaults are 4mm for row and column annotations, and 1cm
#'    for the color legend.
#'    * The `legend_padding` is useful to minimize overlap with
#'    legend and the y-axis labels from the metaplots at the top
#'    of each heatmap.
#' @param profile_value `character` string to define the type of numeric
#'    profile to display at the top of each heatmap. This argument is
#'    passed to `EnrichedHeatmap::anno_enriched()`. Values: `"mean"` the
#'    mean profile; `"sum"` the sum; `"abs_sum"` sum of absolute values;
#'    `"abs_mean"` the mean of absolute values.
#' @param profile_linetype `numeric` or `character` default c(1, 5, 3)
#'    passed to `grid::gpar(lty)`
#'    to define the metaplot line type. Default lty=1 is a solid line.
#'    Values are recycled to the number of profile plots.
#' @param profile_linewidth `numeric` line width passed to `grid::gpar(lwd)`
#'    to control the line width. Default uses lwd=1.5.
#'    Values are recycled to the number of profile plots.
#' @param ylims `numeric` vector of maximum y-axis values for each heatmap
#'    profile; or `list` of min,max values to apply to each `nmatlist` entry.
#' @param border `logical` indicating whether to draw a border around the
#'    heatmap, which includes all heatmap panels in the event of
#'    splitting by clustering. The `border` can be supplied as a vector,
#'    so the `border` can be applied specifically to each heatmap
#'    if needed.
#' @param iter.max `integer` value indicating the maximum iterations
#'    performed by k-means clustering, only relevant when `k_clusters`
#'    is non-zero.
#' @param use_raster `logical` indicating whether to create heatmaps
#'    using raster resizing, almost always recommended `TRUE`
#'    otherwise the output will be very sub-optimal.
#' @param raster_by_magick `logical` passed to `ComplexHeatmap::Heatmap()`,
#'    to enable ImageMagick use during rasterization. By default this
#'    option is `TRUE` and is only disabled when the R package
#'    `"magick"` is not installed, or not properly configured.
#'    If you see a warning "instalilng 'magick' will improve rasterization"
#'    then check the R package with `library(magick)` and see if
#'    there are error messages. When `"magick"` is not available,
#'    the rasterization is substantially slower, and may produce
#'    files much larger than normal.
#' @param raster_quality `logical` passed to `ComplexHeatmap::Heatmap()`,
#'    used when `use_raster=TRUE` and defines the level of detail retained,
#'    and is used only when `raster_by_magick=FALSE`. Using larger numbers
#'    decreases speed substantially.
#' @param do_plot `logical` indicating whether to draw the heatmaps,
#'    `do_plot=TRUE` (default) renders the plots as normal.
#'    `do_plot=FALSE` will return the data used to create heatmaps
#'    without drawing the heatmaps.
#' @param do_caption `logical` indicating whether to include a small caption
#'    at the bottom-right of the plot, describing the number of rows and
#'    columns, the partition, k-means clustering, and main heatmap.
#' @param legend_fontsize `numeric` fontsize to use for all legend text,
#'    default 10.
#'    * Optionally two values can be defined, the first is used for
#'    legend title, the second is used for legend labels.
#' @param padding `grid::unit` object used during `ComplexHeatmap::draw()`
#'    to add whitespace padding around the boundaries of the overall list
#'    of heatmaps. This padding is useful to enforce extra whitespace,
#'    or to prevent labels from exceeding the width of the figure.
#' @param return_type `character` string indicating the type of
#'    data to return:
#'    * `"heatmaplist"` returns the list of heatmaps,
#'    which can separately be arranged together using
#'    `ComplexHeatmap::draw()` or `grid::grid.draw()`.
#'    * `"grid"` returns the `grid` graphical object which may be easier
#'    to render using something like the `patchwork` or `cowplot` R packages.
#' @param show_error `logical` indicating whether to add error
#'    bars to the profile plot at the top of each heatmap.
#'    These error bars are calculated by
#'    `EnrichedHeatmap::anno_enriched()` using
#'    `matrixStats::colSds(x)/nrow(x)`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to
#'    `EnrichedHeatmap::EnrichedHeatmap()` to allow greater
#'    customization of details. Note that many `...` arguments
#'    are also passed to `ComplexHeatmap::Heatmap()`.
#'
#' @family jam coverage heatmap functions
#'
#' @examples
#' ## There is a small example file to use for testing
#' # library(jamba)
#' cov_file1 <- system.file("data", "tss_coverage.matrix", package="platjam");
#' cov_file2 <- system.file("data", "h3k4me1_coverage.matrix", package="platjam");
#' cov_files <- c(cov_file1, cov_file2);
#' names(cov_files) <- gsub("[.]matrix",
#'    "",
#'    basename(cov_files));
#' nmatlist <- coverage_matrix2nmat(cov_files, verbose=FALSE);
#' sapply(nmatlist, function(nmat){attr(nmat, "signal_name")})
#' nmatlist2heatmaps(nmatlist);
#'
#' # sometimes data transform can be helpful
#' nmatlist2heatmaps(nmatlist,
#'    transform=c("log2signed", "sqrt"));
#'
#' # k-means clusters, default uses euclidean distance
#' nmatlist2heatmaps(nmatlist, k_clusters=4,
#'    transform=c("log2signed", "sqrt"));
#'
#' # k-means clusters, "correlation" or "pearson" sometimes works better
#' nmatlist2heatmaps(nmatlist,
#'    k_clusters=4,
#'    min_rows_per_k=20,
#'    k_method="pearson",
#'    transform=c("log2signed", "sqrt"));
#'
#' # example showing usage of top_axis_side
#' # and panel_groups
#' nmatlist2 <- nmatlist[c(1, 1, 1, 2, 2, 2)];
#' names(nmatlist2) <- jamba::makeNames(names(nmatlist2))
#' for (iname in names(nmatlist2)) {
#'    attr(nmatlist2[[iname]], "signal_name") <- gsub("coverage", "cov", iname);
#' }
#' # top_axis_side="left"
#' # assumes 12x7 figure size
#' nmatlist2heatmaps(nmatlist2,
#'    signal_ceiling=0.8,
#'    nmat_colors=rep(c("firebrick", "tomato"), each=3),
#'    panel_groups=rep(c("tss", "h3k4me1"), each=3),
#'    ht_gap=grid::unit(4, "mm"),
#'    top_axis_side="left",
#'    transform=rep(c("log2signed", "sqrt"), each=3));
#'
#' # top_axis_side="both"
#' nmatlist2heatmaps(nmatlist2,
#'    panel_groups=rep(c("tss", "h3k4me1"), each=3),
#'    ht_gap=grid::unit(6, "mm"),
#'    top_axis_side="both",
#'    transform=rep(c("log2signed", "sqrt"), each=3));
#'
#' # multiple heatmap rows
#' nmatlist2heatmaps(nmatlist2,
#'    k_clusters=4,
#'    k_method="pearson",
#'    hm_nrow=2,
#'    panel_groups=rep(c("tss", "h3k4me1"), each=3),
#'    ht_gap=grid::unit(6, "mm"),
#'    top_axis_side="both",
#'    top_anno_height=grid::unit(0.8, "cm"),
#'    transform=rep(c("log2signed", "sqrt"), each=3));
#'
#' # invent anno_df data.frame of additional annotations
#' anno_df <- data.frame(
#'    tss_score=EnrichedHeatmap::enriched_score(jamba::log2signed(nmatlist[[1]])),
#'    h3k4me1_score=EnrichedHeatmap::enriched_score(jamba::log2signed(nmatlist[[2]])),
#'    chromosome=paste0("chr", sample(1:4, replace=TRUE, size=nrow(nmatlist[[1]])))
#' );
#' rownames(anno_df) <- rownames(nmatlist[[1]]);
#' nmatlist2heatmaps(nmatlist,
#'    title="k-means clustering across both heatmaps",
#'    k_clusters=4,
#'    k_method="pearson",
#'    k_heatmap=c(1, 2),
#'    ht_gap=grid::unit(6, "mm"),
#'    top_axis_side="left",
#'    anno_df=anno_df,
#'    transform=rep(c("log2signed", "sqrt"), each=3));
#'
#' # example showing k-means clustering together with annotation groups
#' anno_df <- data.frame(
#'    group=sample(c(1, -1, -1),
#'       size=nrow(nmatlist[[1]]),
#'       replace=TRUE),
#'    row.names=rownames(nmatlist[[1]]))
#' # note for this example the color legends are oriented vertically
#' # showing how the width is adjusted
#' nmatlist2heatmaps(nmatlist,
#'    heatmap_legend_direction="vertical",
#'    k_clusters=0,
#'    color_sub=c(`A`="firebrick", `B`="darkorchid"),
#'    k_colors=c("firebrick", "dodgerblue"),
#'    min_rows_per_k=50,
#'    ht_gap=grid::unit(1, "cm"),
#'    k_method="correlation",
#'    k_heatmap=1:2,
#'    anno_df=anno_df,
#'    partition="group",
#'    row_title_rot=0,
#'    transform=rep(c("log2signed", "sqrt"), each=3));
#'
#' # same as above, partition and k_clusters together
#' # except uses multiple values for k_clusters
#' nmatlist2heatmaps(nmatlist,
#'    k_clusters=c(1, 4),
#'    min_rows_per_k=25,
#'    k_heatmap=1:2,
#'    k_method="correlation",
#'    anno_df=anno_df,
#'    partition="group",
#'    row_title_rot=0)
#'
#' @export
nmatlist2heatmaps <- function
(nmatlist,
 panel_groups=NULL,
 title=NULL,
 title_gp=grid::gpar(fontsize=16),
 caption=NULL,
 upstream_length=NULL,
 downstream_length=NULL,
 k_clusters=0,
 min_rows_per_k=100,
 k_subset=NULL,
 k_colors=NULL,
 k_width=grid::unit(5, "mm"),
 k_method=c("correlation",
    "euclidean",
    "pearson",
    "spearman"),
 k_heatmap=main_heatmap,
 partition=NULL,
 row_title_rot=0,
 partition_counts=TRUE,
 partition_count_template="{partition_name}\n({counts} rows)",
 rows=NULL,
 row_order=NULL,
 nmat_colors=NULL,
 middle_color="white",
 nmat_names=NULL,
 main_heatmap=NULL,
 anno_df=NULL,
 byCols=NULL,
 color_sub=NULL,
 anno_row_marks=NULL,
 anno_row_labels=NULL,
 anno_row_gp=grid::gpar(fontsize=14),
 recenter_heatmap=NULL,
 summit_names=NULL,
 recenter_range=NULL,
 recenter_invert=FALSE,
 restrand_heatmap=NULL,
 restrand_range=NULL,
 restrand_buffer=NULL,
 restrand_invert=FALSE,
 top_annotation=NULL,
 top_anno_height=grid::unit(3, "cm"),
 top_axis_side=c("right"),
 legend_max_ncol=2,
 legend_base_nrow=12,
 legend_max_labels=40,
 show_heatmap_legend=TRUE,
 heatmap_legend_param=NULL,
 heatmap_legend_direction="horizontal",
 annotation_legend_param=NULL,
 hm_nrow=1,
 transform="none",
 transform_label=NULL,
 signal_ceiling=NULL,
 axis_name=NULL,
 axis_name_gp=grid::gpar(fontsize=10),
 axis_name_rot=90,
 column_title_gp=grid::gpar(fontsize=14),
 lens=-2,
 anno_lens=8,
 pos_line=FALSE,
 seed=123,
 ht_gap=grid::unit(4, "mm"),
 row_anno_padding=grid::unit(4, "mm"),
 column_anno_padding=grid::unit(4, "mm"),
 legend_padding=grid::unit(1, "cm"),
 profile_value=c("mean",
    "sum",
    "abs_mean",
    "abs_sum"),
 profile_linetype=c(1, 5, 3),
 profile_linewidth=1.5,
 ylims=NULL,
 border=TRUE,
 iter.max=20,
 use_raster=TRUE,
 raster_quality=1,
 raster_by_magick=jamba::check_pkg_installed("magick"),
 do_plot=TRUE,
 do_caption=TRUE,
 legend_fontsize=10,
 legend_width=grid::unit(3, "cm"),
 trim_legend_title=TRUE,
 padding=grid::unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
 return_type=c("heatmaplist", "grid"),
 show_error=FALSE,
 verbose=FALSE,
 ...)
{
   #
   return_type <- match.arg(return_type);
   profile_value <- match.arg(profile_value);
   if (length(legend_fontsize) == 0) {
      legend_fontsize <- 10;
   }
   legend_fontsize <- rep(legend_fontsize, length.out=2);
   if (length(seed) > 0) {
      set.seed(seed);
   }
   if (length(main_heatmap) == 0) {
      main_heatmap <- seq_along(nmatlist)
   }
   if (any(main_heatmap > length(nmatlist))) {
      main_heatmap <- main_heatmap[main_heatmap > 0 &
            !is.na(main_heatmap) &
            main_heatmap %in% seq_along(nmatlist)];
   }

   if (length(border) == 0) {
      border <- FALSE;
   }
   border <- rep(border, length.out=length(nmatlist));
   if (length(legend_width) == 0) {
      legend_width <- grid::unit(3, "cm");
   }

   ##############################################################
   # Optionally re-center rows
   if (length(recenter_heatmap) > 0 || length(summit_names) > 0) {
      nmatlist <- recenter_nmatlist(nmatlist=nmatlist,
         recenter_heatmap=recenter_heatmap,
         summit_names=summit_names,
         recenter_range=recenter_range,
         recenter_invert=recenter_invert,
         verbose=verbose,
         ...);
   }

   ##############################################################
   # Optionally re-strand rows
   if (length(restrand_heatmap) > 0) {
      nmatlist <- restrand_nmatlist(nmatlist=nmatlist,
         restrand_heatmap=restrand_heatmap,
         restrand_range=restrand_range,
         restrand_buffer=restrand_buffer,
         restrand_invert=restrand_invert,
         verbose=verbose,
         ...);
   }
   # Summarize recenter/restrand adjustments
   adjust_df <- NULL;
   if (length(recenter_heatmap) > 0 || length(restrand_heatmap) > 0) {
      nmatlist_dfs <- nmatlist_summary(nmatlist);
      adjust_df <- jamba::mergeAllXY(nmatlist_dfs);
      # print(head(adjust_df, 30));# debug
   }

   ##############################################################
   # optional coordinate range zoom for coverage data
   if (length(upstream_length) > 0 || length(downstream_length) > 0) {
      if (verbose > 1) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "call zoom_nmatlist()");
      }
      nmatlist <- zoom_nmatlist(nmatlist=nmatlist,
         upstream_length=upstream_length,
         downstream_length=downstream_length);
   }

   # store original (default) heatmap annotation padding
   ROW_ANNO_PADDING <- ComplexHeatmap::ht_opt("ROW_ANNO_PADDING")
   COLUMN_ANNO_PADDING <- ComplexHeatmap::ht_opt("COLUMN_ANNO_PADDING")
   HEATMAP_LEGEND_PADDING <- ComplexHeatmap::ht_opt("HEATMAP_LEGEND_PADDING")
   ANNOTATION_LEGEND_PADDING <- ComplexHeatmap::ht_opt("ANNOTATION_LEGEND_PADDING")

   # define custom row/column heatmap annotation padding
   if (length(row_anno_padding) == 0) {
      row_anno_padding <- grid::unit(4, "mm");
   }
   if (length(column_anno_padding) == 0) {
      column_anno_padding <- grid::unit(4, "mm");
   }
   if (!grid::is.unit(row_anno_padding)) {
      row_anno_padding <- grid::unit(row_anno_padding, "mm")
   }
   if (!grid::is.unit(column_anno_padding)) {
      column_anno_padding <- grid::unit(column_anno_padding, "mm")
   }
   # define custom row/column heatmap annotation padding
   ComplexHeatmap::ht_opt("ROW_ANNO_PADDING"=row_anno_padding)
   ComplexHeatmap::ht_opt("COLUMN_ANNO_PADDING"=column_anno_padding)
   ComplexHeatmap::ht_opt(HEATMAP_LEGEND_PADDING=legend_padding)
   ComplexHeatmap::ht_opt(ANNOTATION_LEGEND_PADDING=legend_padding)
   ## k_method
   kmeans <- function(centers, ...){stats::kmeans(centers=centers, ...)};
   # k_method <- head(k_method, 1);
   k_method <- match.arg(k_method);
   if (!jamba::check_pkg_installed("amap")) {
      k_method <- "euclidean";
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "k_method requires the '",
            "amap",
            "' package, which is not installed. Setting ",
            "k_method='euclidean'");
      }
   }
   if (!"euclidean" %in%  k_method) {
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "Using amap::Kmeans()");
      }
      kmeans <- function(centers, ...){
         amap::Kmeans(...,
            centers=centers,
            method=k_method)
      };
   }
   nmat_rows <- Reduce("intersect",
      lapply(nmatlist, rownames));
   if (length(nmat_rows) == 0) {
      stop("There were no shared rownames across data matrices in nmatlist.");
   }
   if (length(rows) == 0) {
      ## Make sure rows are present in all nmatlist entries.
      rows <- nmat_rows;
   } else {
      ## Make sure rows are present in all rownames of nmatlist
      if (is.numeric(rows)) {
         rows <- rmNA(nmat_rows[rows]);
      } else {
         rows <- intersect(rows, nmat_rows);
         # rows <- rows[rows %in% nmat_rows];
      }
      if (length(rows) == 0) {
         stop("No values in rows matched any nmatlist rownames.");
      }
   }
   ## Also optionally subset rows by rownames(anno_df)
   if (length(anno_df) > 0) {
      rows <- intersect(rows, rownames(anno_df));
      # rows <- rows[rows %in% rownames(anno_df)];
      if (length(rows) == 0) {
         stop("rownames(anno_df) did not match any values in rows.");
      }
   }
   partition_colors <- NULL;
   if (length(partition) > 0) {
      if (length(anno_df) > 0 && all(partition %in% colnames(anno_df))) {
         partition_df <- data.frame(check.names=FALSE,
            anno_df[, partition, drop=FALSE]);
         for (pcol in partition) {
            if (is.numeric(partition_df[[pcol]])) {
               # sort numeric columns decreasing order
               # so numbers will be ordered bottom to top (like scatterplots)
               partition_df[[pcol]] <- factor(partition_df[[pcol]],
                  levels=rev(sort(unique(partition_df[[pcol]]))))
            }
         }
         partition <- jamba::pasteByRowOrdered(partition_df)
         names(partition) <- rownames(anno_df);
      }
      if (length(names(partition)) > 0) {
         rows <- intersect(rows, names(partition));
         # rows <- rows[rows %in% names(partition)];
         if (length(rows) == 0) {
            stop("names(partition) did not match any rownames in nmatlist.");
         }
         partition <- partition[match(rows, names(partition))];
      } else {
         stop("names(partition) were not provided, and must match rownames in nmatlist.");
      }
      # try to intuit partition colors
      if (!is.factor(partition)) {
         partition <- factor(partition,
            levels=jamba::mixedSort(unique(partition)))
      } else {
         # call factor() which removes empty factor levels
         partition <- factor(partition);
      }
      partition_values <- levels(factor(partition));
      if (length(color_sub) > 0) {
         if (is.atomic(color_sub)) {
            if (all(partition_values %in% names(color_sub))) {
               partition_colors <- color_sub[partition_values];
            }
         } else if ("color_sub" %in% names(attributes(color_sub))) {
            color_sub_v <- attr(color_sub, "color_sub");
            if (all(partition_values %in% names(color_sub_v))) {
               partition_colors <- color_sub_v[partition_values];
            }
         }
         if (length(partition_colors) == 0 && is.atomic(color_sub)) {
            if (all(partition_values %in% names(color_sub))) {
               partition_colors <- color_sub[partition_values];
            }
         }
         if (length(partition_colors) == 0 && is.list(color_sub)) {
            for (icol in names(color_sub)) {
               if (length(partition_colors) == 0 &&
                     is.atomic(names(color_sub[[icol]])) &&
                     all(partition_values %in% names(color_sub[[icol]]))) {
                  partition_colors <- color_sub[[icol]][partition_values];
               }
            }
         }
      }
      # fallback to categorical assignment
      if (length(partition_colors) == 0) {
         if (length(k_colors) >= length(partition_values)) {
            partition_colors <- head(k_colors, length(partition_values));
            names(partition_colors) <- partition_values;
         } else {
            partition_colors <- colorjam::group2colors(x=partition_values,
               sortFunc=c)
         }
      }
      # replace any NA colors with grey?
      if (any(is.na(partition_colors))) {
         partition_colors[is.na(partition_colors)] <- "#AAAAAA"
      }
   }

   #########################################
   ## Summarize rows recognized thus far
   if (verbose) {
      jamba::printDebug("nmatlist2heatmaps(): ",
         "Recognized ",
         jamba::formatInt(length(rows)),
         " rows shared across all matrices.");
   }

   #########################################
   ## Optional transform for each matrix
   if (length(transform) == 0) {
      transform <- function(x){x}
   }
   if (length(transform_label) == 0) {
      transform_label <- sapply(seq_along(transform), function(itr1){
         if (length(names(transform)) > 0 && all(nchar(names(transform)[itr1]) > 0)) {
            names(transform)[itr1];
         } else if (is.character(transform[[itr1]])){
            transform[[itr1]]
         } else {
            ""
         }
      })
   }
   transform_label <- rep(transform_label, length.out=length(nmatlist))
   # convert transform to proper function
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

   #########################################
   ## validate row_order
   row_order_type <- NULL;
   if (length(row_order) == 0) {
      row_order <- TRUE;
   }
   ## row_order logical TRUE or FALSE
   if (is.logical(row_order)) {
      if (TRUE %in% row_order) {
         row_order_type <- "Enrichment Score";
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ", sep="",
               c("Define row_order for main_heatmap (",
                  jamba::cPasteS(main_heatmap),
                  ") using ",
                  "EnrichedHeatmap::enriched_score()"));
         }
         ## Calculate scores across each main_heatmap
         ## Todo: Consider applying data transform prior to this step.
         row_scores <- lapply(main_heatmap, function(i){
            row_score <- jamba::rmNA(
               naValue=0,
               EnrichedHeatmap::enriched_score(
                  transform[[i]](
                     nmatlist[[i]][rows, , drop=FALSE])
               )
            )
         })
         row_score <- jamba::rmNA(naValue=0,
            Reduce("+", row_scores));
         row_order <- order(row_score, decreasing=TRUE);
         names(row_order) <- rows;
      } else {
         row_order <- FALSE;
      }
   } else if (length(row_order) > 1) {
      ## row_order must be named by rows
      if (length(row_order) < length(rows)) {
         cli::cli_abort(message=paste0(
            "To few {.var row_order} values were supplied",
            " ({length(row_order)}),",
            " for the number of rownames ({length(rows)})."))
      }
      if (!is.numeric(row_order)) {
         cli::cli_abort(message=paste0(
            "{.var row_order} must be supplied as a {.cls numeric} vector,",
            " preferably using {.code names(row_order)} to match rownames."))
      }
      if (any(duplicated(row_order))) {
         cli::cli_abort(message=paste0(
            "{.var row_order} contained duplicate values which are not valid.",
            " If these values are actually scores, they can be converted to",
            " order using {.code order(row_order)}."));
      }
      if (!all(row_order == round(row_order))) {
         cli::cli_abort(message=paste0(
            "{.var row_order} contained non-integer values which are not",
            " valid. If these values are actually scores, they can be ",
            " converted toorder using {.code order(row_order)}."));
      }
      if (length(names(row_order)) == 0 &&
            length(row_order) == length(rows)) {
         cli::cli_warn(message=paste0(
            "{.var row_order} was supplied as {.cls numeric} values without",
            " names.",
            " Supply {.code names(row_order)} to guarantee ",
            "correct ordering."))
         names(row_order) <- rows;
      } else if (!all(rows %in% names(row_order))) {
         cli::cli_abort(message=paste0(
            "{.var row_order} supplied as {.cls numeric} values did not",
            " not represent all rownames in {.code names(row_order)}."))
      } else {
         row_order <- row_order[match(rows, names(row_order))];
         row_order_type <- "Provided Row Order";
      }
   }

   #########################################
   ## validate panel_groups
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

   if (FALSE) {
      ## Optional transformation of each matrix
      if (length(transform) == 0) {
         transform <- function(x){x}
      }
      if (length(transform_label) == 0) {
         transform_label <- sapply(seq_along(transform), function(itr1){
            if (length(names(transform)) > 0 && all(nchar(names(transform)[itr1]) > 0)) {
               names(transform)[itr1];
            } else if (is.character(transform[[itr1]])){
               transform[[itr1]]
            } else {
               ""
            }
         })
      }
      transform_label <- rep(transform_label, length.out=length(nmatlist))
      # convert transform to proper function
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
   }

   ## pos_line
   if (length(pos_line) == 0) {
      pos_line <- FALSE;
   }
   pos_line <- rep(pos_line,
      length.out=length(nmatlist));

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

   ##################################
   ## Optional k-means clustering
   if (length(k_clusters) > 0 && any(k_clusters > 1)) {
      # expand k_clusters to the number of partitions
      if (length(partition) > 0) {
         k_clusters <- rep(k_clusters,
            length.out=length(unique(partition)))
      } else {
         k_clusters <- head(k_clusters, 1)
      }
      # 0.0.76.900 - enforce maximum k_clusters value using min_rows_per_l
      k_clusters_max <- ceiling(length(rows) / min_rows_per_k);
      k_clusters <- jamba::noiseFloor(k_clusters,
         minimum=1,
         ceiling=min(k_clusters_max))

      if (length(k_clusters) == 1) {
         # Note: only assign colors if there are no partitions
         if (length(k_colors) == 0) {
            k_colors <- colorjam::group2colors(
               seq_len(k_clusters),
               colorSub=color_sub);
         } else if (length(k_colors) < k_clusters) {
            ## Expand the given colors using color2gradient()
            k_multiplier <- ceiling(k_clusters / length(k_colors));
            k_colors <- jamba::nameVector(
               rev(head(
                  jamba::color2gradient(k_colors,
                     n=k_multiplier,
                     gradientWtFactor=2/3),
                  k_clusters)),
               seq_len(k_clusters));
         } else if (length(names(k_colors)) == 0) {
            names(k_colors) <- rev(seq_len(k_clusters));
         }
      }
      if (length(k_heatmap) == 0) {
         k_heatmap <- main_heatmap;
      }
      # for multiple k_heatmap values, cbind each transformed matrix
      if (length(k_heatmap) > 1) {
         imatrix <- do.call(cbind, lapply(k_heatmap, function(k_heatmap1){
            kmatch <- match(rows, rownames(nmatlist[[k_heatmap1]]))
            transform[[k_heatmap1]](
               nmatlist[[k_heatmap1]][kmatch, , drop=FALSE]);
         }));
      } else {
         kmatch <- match(rows, rownames(nmatlist[[k_heatmap]]))
         imatrix <- transform[[k_heatmap]](
            nmatlist[[k_heatmap]][kmatch, , drop=FALSE]);
      }
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "Running kmeans, k_clusters:",
            k_clusters,
            ", k_method:",
            k_method);
      }

      # kmeans cluster within each partition when defined
      if (length(partition) == 0) {
         kpartition <- suppressMessages(suppressWarnings(
            kmeans(
               imatrix,
               iter.max=iter.max,
               centers=head(k_clusters, 1))$cluster));
         ## Confirm that names(partition) match rows
         names(kpartition) <- rows;
         partition <- kpartition;
         if (verbose) {
            k_sizes <- table(kpartition);
            jamba::printDebug("nmatlist2heatmaps(): ",
               "k-means cluster sizes: ",
               paste0("cluster", names(k_sizes), "=", k_sizes), sep=", ");
         }
      } else {
         ## previous strategy for row partitions and k-means clusters:
         # - perform global k-means, then split by partition then k-cluster
         # - downside is each partition is not individually clustered
         #   so any one k-cluster may only have 1 or 2 members in the partition
         ## current strategy:
         # - iterate each partition, apply k-means within the partition
         # - each partition can thus have different k values, based upon size
         # - or k values recycled from user-provided k_clusters
         partition_rows_list <- split(rows, partition[rows]);
         # recycle k_clusters across each partition
         # if k_clusters has names that match partition_rows_list,
         # allow name-directed assignment
         if (all(names(partition_rows_list) %in% names(k_clusters))) {
            k_clusters <- k_clusters[names(partition_rows_list)];
         } else {
            k_clusters <- rep(k_clusters,
               length.out=length(partition_rows_list));
            names(k_clusters) <- names(partition_rows_list);
         }
         # requires min_rows_per_k rows per cluster by default, or adjust k
         kpartitions_list <- lapply(seq_along(partition_rows_list), function(ipart){
            prows <- partition_rows_list[[ipart]];
            use_k <- k_clusters[[ipart]];
            if (use_k > ceiling(length(prows) / min_rows_per_k)) {
               use_k <- ceiling(length(prows) / min_rows_per_k);
            }
            if (length(prows) <= min_rows_per_k || use_k %in% c(0, 1, NA)) {
               return(jamba::nameVector(
                  rep("", length.out=length(prows)),
                  prows))
            }
            kpartition <- suppressMessages(suppressWarnings(
               kmeans(
                  imatrix[match(prows, rownames(imatrix)), , drop=FALSE],
                  iter.max=iter.max,
                  centers=use_k)$cluster));
            names(kpartition) <- prows;
            kpartition;
         })
         kpartition <- jamba::rmNA(naValue=0,
            unlist(unname(kpartitions_list))[rows]);
         names(kpartition) <- rows;
         partition_df <- data.frame(
            partition=partition[rows],
            kpartition=kpartition[rows]);
         # create new partition that honors factor order where relevant
         partition_sep <- ": ";
         partition <- jamba::nameVector(
            jamba::pasteByRowOrdered(partition_df,
               sep=partition_sep),
            rows);
         # assign partition_colors using color2gradient to split colors
         # by each partition
         partition_dfu <- jamba::mixedSortDF(unique(partition_df));
         use_partition_colors <- partition_colors[
            as.character(partition_dfu$partition)];
         if (any(duplicated(use_partition_colors))) {
            use_partition_colors <- jamba::color2gradient(use_partition_colors,
               dex=2)
         }
         names(use_partition_colors) <- jamba::pasteByRow(partition_dfu,
            sep=partition_sep);
         partition_colors <- use_partition_colors;
         k_colors <- partition_colors;
         # k_colors <- colorjam::group2colors(levels(partition),
         #    colorSub=color_sub);
      }
   }
   if (verbose > 1 && length(k_colors) > 0) {
      jamba::printDebug("post k-means k_colors: ");jamba::printDebugI(k_colors);# debug
   }
   # End k-means clustering
   ##################################

   ##################################
   ## Partition heatmap sidebar
   if (length(partition) > 0) {
      ## Make sure to use the partition values with the properly ordered rows
      if (!all(rows %in% names(partition))) {
         rows <- intersect(rows, names(partition));
         # rows <- rows[rows %in% names(partition)];
      }
      if (!any(rows %in% names(partition))) {
         missing_rows <- setdiff(rows, names(partition));
         jamba::printDebug("missing_rows (", jamba::formatInt(length(missing_rows)), "):");
         print(head(missing_rows));
         jamba::printDebug("head(partition):");
         print(head(partition));
         jamba::printDebug("head(rows):");
         print(head(rows));
         print(table(all(rows %in% names(partition))));
         stop("names(partition) must match rownames in nmatlist.");
      }
      partition <- partition[match(rows, names(partition))];
      if (!is.factor(partition)) {
         partition <- factor(partition,
            levels=jamba::mixedSort(unique(partition)));
      } else {
         # remove empty factor levels
         partition <- factor(partition);
      }
      ## Define colors if not provided
      if (length(k_colors) == 0) {
         k_colors <- jamba::nameVector(
            # colorjam::rainbowJam(length(levels(partition)), ...),
            partition_colors[levels(partition)],
            levels(partition));
         # k_colors <- k_colors[sort(names(k_colors))];
      } else {
         if (length(names(k_colors)) == 0) {
            if (length(k_colors) >= length(levels(partition))) {
               k_colors <- head(k_colors, length(levels(partition)));
               names(k_colors) <- levels(partition);
            } else {
               k_colors <- colorjam::rainbowJam(length(levels(partition)), ...)
               names(k_colors) <- levels(partition);
            }
         }
         if (all(levels(partition) %in% names(k_colors))) {
            k_colors <- k_colors[levels(partition)]
         }
      }
      if (verbose > 1 && length(k_colors) > 0) {
         jamba::printDebug("partition step 1 k_colors: ");jamba::printDebugI(k_colors);# debug
      }

      ## Optional subset of k-means clusters
      if (length(k_subset) > 0) {
         if (!any(as.character(partition) %in% as.character(k_subset))) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Warning: k_subset was supplied but does not match any partition values.");
            jamba::printDebug("nmatlist2heatmaps(): ",
               "head(levels(partition), 20):",
               paste0("'", head(levels(partition), 20), "'"));
         }
         partition_keep <- (as.character(partition) %in% as.character(k_subset));
         partition <- factor(partition[partition_keep],
            levels=intersect(as.character(k_subset),
               levels(partition)));
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
      p_num <- length(levels(partition[rows]));
      p_ncol <- min(c(ceiling(p_num / legend_base_nrow), legend_max_ncol));
      p_nrow <- ceiling(p_num / p_ncol);
      # p_at <- jamba::mixedSort(unique(partition[rows]));
      p_at <- levels(factor(partition[rows]))
      ## p_labels is defined without row counts
      # p_labels <- gsub("\n", " ", p_at);
      p_labels <- p_at;

      ## Optionally rename partitions to include row counts
      # partition_counts <- TRUE;
      # partition_count_template <- "{partition_name}\n({counts} rows)";
      if (TRUE %in% partition_counts) {
         partition_cts <- table(factor(partition[rows]));
         p_after <- sapply(seq_along(partition_cts), function(ip){
            partition_name <- p_at[ip];
            counts <- partition_cts[[partition_name]];
            glue::glue(partition_count_template,
               partition_name=partition_name,
               counts=jamba::formatInt(counts))
         })
         p_before_after <- data.frame(before=p_at,
            after=p_after)
         level_match <- match(p_before_after$before, levels(partition))
         levels(partition)[level_match] <- p_before_after$after;
         level_k_match <- match(p_before_after$before, names(k_colors))
         names(k_colors)[level_k_match] <- p_before_after$after;
         p_at <- p_before_after$after;
      }

      p_heatmap_legend_param <- list(
         title="Cluster",
         title_position="topleft",
         type="boxplot",
         pch=26:28,
         border="black",
         nrow=p_nrow,
         at=p_at,
         labels=p_labels,
         labels_gp=grid::gpar(
            fontsize=legend_fontsize[2]),
         title_gp=grid::gpar(
            fontsize=legend_fontsize[1],
            fontface="bold"),
         legend_gp=grid::gpar(
            lty=profile_linetype,
            lwd=profile_linewidth)
      )
      PHM <- ComplexHeatmap::Heatmap(partition[rows],
         border=FALSE,
         heatmap_legend_param=p_heatmap_legend_param,
         use_raster=use_raster,
         column_names_rot=axis_name_rot,
         raster_quality=raster_quality,
         raster_by_magick=raster_by_magick,
         col=k_colors,
         name="cluster",
         show_row_names=FALSE,
         width=k_width);
      PHM_rows <- rows;
   }


   ##########################################################
   ## Optional data.frame with additional annotations
   if (length(anno_df) > 0) {
      if (!jamba::igrepHas("data.frame|dataframe|data.table|tibble", class(anno_df))) {
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
         if (is.numeric(row_order) && !any(duplicated(row_order))) {
            row_rank <- match(rows, rows[row_order]);
            anno_df <- jamba::mixedSortDF(
               data.frame(check.names=FALSE,
                  anno_df[jamba::rmNA(match(rows, rownames(anno_df))),, drop=FALSE],
                  row_rank_JAM=row_rank),
               byCols=c(byCols, "row_rank_JAM"));
            anno_df <- anno_df[,setdiff(
               colnames(anno_df),
               "row_rank_JAM"), drop=FALSE];
            byCols <- c(byCols,
               row_order_type);
         } else {
            anno_df <- jamba::mixedSortDF(anno_df,
               byCols=byCols);
         }
         rows <- intersect(rownames(anno_df), rows);
         # rows <- rownames(anno_df)[rownames(anno_df) %in% rows];
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Sorted rows by:",
               byCols);
         }
         row_order <- FALSE;
      } else {
         rows <- rows[rows %in% rownames(anno_df)];
      }
      anno_df <- anno_df[match(rows, rownames(anno_df)),, drop=FALSE];
      ## Determine a list of color functions, one for each column
      anno_colors_l <- lapply(jamba::nameVector(colnames(anno_df)), function(i){
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "anno_colors_l colname:", i);
         }
         i1 <- jamba::rmNA(anno_df[[i]]);
         if (any(c("integer", "numeric") %in% class(i1))) {
            if (min(i1, na.rm=TRUE) < 0 || max(abs(i1), na.rm=TRUE) < 50) {
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
                  lens=anno_lens,
                  trimRamp=c(2, 2),
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
               iminmax <- quantile(i1,
                  c(0.005, 0.995));
               ibreaks <- unique(seq(from=iminmax[1],
                  to=iminmax[2],
                  length.out=15));
               if (max(ibreaks) == 0 || length(ibreaks) <= 1) {
                  ibreaks <- c(0, 1);
               }
               colBR <- jamba::getColorRamp("Purples",
                  n=length(ibreaks),
                  lens=anno_lens);
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
            i2 <- jamba::mixedSort(unique(jamba::rmNA(i1)));
            if (all(jamba::isColor(i2))) {
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "anno_colors_l colname:", i,
                     " pre-defined colors");
               }
               cBR <- jamba::nameVector(i2);
            } else if (length(color_sub) > 0 && all(i2 %in% names(color_sub))) {
               color_match <- match(i2, names(color_sub));
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "anno_colors_l colname:", i,
                     " categorical data using color_sub");
               }
               cBR <- color_sub[color_match];
            } else {
               if (!"factor" %in% class(i2)) {
                  i2 <- factor(i2,
                     levels=jamba::mixedSort(i2));
               }
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "anno_colors_l colname:", i,
                     " categorical data");
               }
               cBR <- colorjam::group2colors(
                  i2,
                  colorSub=color_sub);
               cBR <- jamba::rmNA(cBR);
            }
         }
         cBR;
      });
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "sdim(anno_colors_l):");
         print(jamba::sdim(anno_colors_l));
         # str(anno_colors_l);
         platjam::print_color_list(anno_colors_l);
      }

      ## annotation_legend_param
      if (length(annotation_legend_param) == 0) {
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Defining annotation_legend_param.");
         }
         ## list show_legend
         annotation_show_legend <- sapply(anno_colors_l, function(ac1){
            if (!is.function(ac1) & length(ac1) > legend_max_labels) {
               FALSE;
            } else {
               TRUE;
            }
         });
         ## list of annotation_legend_param
         annotation_legend_param <- lapply(jamba::nameVector(colnames(anno_df)), function(i){
            i1 <- jamba::rmNA(anno_df[[i]]);
            a_num <- length(setdiff(unique(i1), c("", NA)));
            a_ncol <- min(c(ceiling(a_num / legend_base_nrow), legend_max_ncol));
            a_nrow <- ceiling(a_num / a_ncol);
            i_title <- apply_word_wrap(i, width=15, sep="\n")
            if (a_num <= 10) {
               ## display distinct steps
               if (is.numeric(i1)) {
                  i1_at <- sort(unique(i1));
                  i1_labels <- jamba::formatInt(i1_at);
               } else {
                  i1_at <- jamba::mixedSort(unique(i1));
                  i1_labels <- gsub("\n", " ", i1_at);
               }
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "annotation_legend_param colname:", i,
                     " discrete numeric color legend, a_num:",
                     a_num,
                     ", a_nrow:",
                     a_nrow);
               }
               list(
                  title=i_title,
                  title_position="topleft",
                  at=i1_at,
                  labels=i1_labels,
                  color_bar="discrete",
                  border="black",
                  nrow=a_nrow
               );
            } else if (any(c("integer", "numeric") %in% class(i1))) {
               ## display continuous color gradient
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "annotation_legend_param colname:", i,
                     " continuous color legend");
               }
               list(
                  direction="horizontal",
                  title=i_title,
                  labels_rot=90,
                  legend_width=legend_width,
                  title_position="topleft",
                  border="black",
                  grid_width=grid::unit(1, "npc"));
            } else {
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "annotation_legend_param colname:", i,
                     " discrete categorical color legend");
               }
               if (is.numeric(i1)) {
                  i1_at <- sort(unique(i1));
                  i1_labels <- jamba::formatInt(i1_at);
               } else {
                  i1_at <- jamba::mixedSort(unique(i1));
                  i1_labels <- gsub("\n", " ", i1_at);
               }
               list(
                  title=i_title,
                  title_position="topleft",
                  at=i1_at,
                  labels=i1_labels,
                  border="black",
                  nrow=a_nrow
               )
               #   grid_width=grid::unit(1, "npc")
            }
         });
      }
      for (jj in 1:ncol(anno_df)) {
         i1 <- anno_df[[jj]];
         if (any(c("character") %in% class(i1))) {
            i1 <- factor(i1,
               levels=jamba::mixedSort(unique(i1)));
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
         df=anno_df[rows, , drop=FALSE],
         annotation_legend_param=annotation_legend_param,
         show_legend=annotation_show_legend,
         name="Annotation",
         annotation_name_rot=axis_name_rot,
         col=anno_colors_l);
      AHM_rows <- rows;

      ## Optional row marks
      anno_rows <- rows[rows %in% anno_row_marks];
      if (length(anno_rows) > 0) {
         anno_row_which <- match(anno_rows, rows);
         if (length(anno_row_labels) > 0 && all(anno_row_labels %in% colnames(anno_df))) {
            anno_row_labels <- pasteByRow(
               anno_df[anno_rows,anno_row_labels,drop=FALSE],
               sep=" ");
         } else if (length(anno_row_labels) >= length(anno_rows)) {
            if (length(names(anno_row_labels)) == 0 &&
                  length(anno_row_labels) == length(anno_row_marks)) {
               names(anno_row_labels) <- anno_row_marks;
            }
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
         if (!"gpar" %in% class(anno_row_gp)) {
            anno_row_gp <- grid::gpar(fontsize=14)
         }
         MHM <- ComplexHeatmap::Heatmap(jamba::nameVector(anno_df[rows,1], rows),
            col=anno_colors_l[[1]],
            name=colnames(anno_df)[1],
            show_row_names=FALSE,
            width=k_width,
            cluster_rows=FALSE,
            right_annotation=ComplexHeatmap::rowAnnotation(
               foo=ComplexHeatmap::anno_mark(at=anno_row_which,
                  labels_gp=anno_row_gp,
                  labels=anno_row_labels)
            )
         );
         MHM_rows <- rows;
      } else {
         MHM <- NULL;
      }
   } else {
      AHM <- NULL;
      MHM <- NULL;
      byCols <- row_order_type;
   }


   ##################################
   ## panel_groups
   if (length(panel_groups) > 0) {
      ## Make sure we have some duplicated panel_groups
      if (length(jamba::tcount(panel_groups, minCount=2)) > 0) {
         if (length(show_heatmap_legend) <= 1) {
            if (length(show_heatmap_legend) == 0) {
               show_heatmap_legend <- TRUE;
            }
            show_heatmap_legend <- ifelse(duplicated(panel_groups),
               FALSE,
               show_heatmap_legend);
         }
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
         # Advanced todo: Consider calculating the std error so
         # the y-axis can be expanded to accomodate added range
         # when show_error=TRUE.
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
                  ## jamba::rmNULL() removes empty elements if factor levels
                  ## are not present in the rows of data being used
                  plist <- jamba::rmNULL(split(names(partition),
                     partition));
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
         if (length(ylims) == 0) {
            ylims <- panel_ylims[panel_groups];
         }
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
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "panel_ceilings[panel_groups]: ",
               paste0("(", jamba::cPaste(panel_ceilings[panel_groups]), ")"),
               sep="; ");
         }
      }
   }
   if (length(show_heatmap_legend) == 0) {
      show_heatmap_legend <- TRUE;
   }
   if (length(show_heatmap_legend) != length(nmatlist)) {
      show_heatmap_legend <- rep(show_heatmap_legend,
         length.out=length(nmatlist));
   }

   #############################################
   ## Iterate each matrix to create heatmaps
   if (length(lens) == 0) {
      lens <- 0;
   }
   lens <- rep(lens,
      length.out=length(nmatlist));
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
   if (length(axis_name) == 0) {
      axis_name <- lapply(nmatlist, function(nmat){
         extend1 <- signif(attr(nmat, "extend")[1], digits=2);
         extend2 <- signif(attr(nmat, "extend")[2], digits=2);
         # off-by-one is assumed to be identical
         if (extend1 == (extend2 + 1)) {
            extend2 <- extend2 + 1;
         }
         c(
            paste0("-", jamba::formatInt(extend1)),
            attr(nmat, "target_name"),
            jamba::formatInt(extend2))
      })
   } else if (is.list(axis_name)) {
      axis_name <- rep(axis_name,
         length.out=length(nmatlist));
   } else {
      if (length(axis_name) == 3) {
         axis_name <- rep(list(axis_name),
            length.out=length(nmatlist));
      } else {
         axis_name <- rep(axis_name,
            length.out=length(nmatlist));
         axis_name <- lapply(seq_along(nmatlist), function(inmat){
            nmat <- nmatlist[[inmat]];
            c(
               paste0("-",
                  jamba::formatInt(attr(nmat, "extend")[1])),
               axis_name[[inmat]],
               jamba::formatInt(attr(nmat, "extend")[2]))
         })
      }
   }


   ###############################################
   # expand heatmap_legend_param to each heatmap
   if (length(heatmap_legend_param) == 0) {
      # prototype default heatmap legend
      # heatmap_legend_direction <- "horizontal";
      if ("horizontal" %in% heatmap_legend_direction) {
         hml_grid_width <- grid::unit(1, "npc")
      } else {
         hml_grid_width <- grid::unit(5, "mm")
      }
      heatmap_legend_param_1 <- list(
         direction=heatmap_legend_direction,
         legend_width=legend_width,
         title_position="topleft",
         labels_gp=grid::gpar(
            fontsize=legend_fontsize[2]),
         title_gp=grid::gpar(
            fontsize=legend_fontsize[1],
            fontface="bold"),
         border="black",
         grid_width=hml_grid_width);

      heatmap_legend_param <- lapply(seq_along(nmatlist), function(ipanel){
         hlp_1 <- heatmap_legend_param_1;
         if (length(panel_groups) > 0) {
            hlp_1$title <- panel_groups[[ipanel]];
         } else {
            hlp_1$title <- attr(nmatlist[[ipanel]], "signal_name");
         }
         if (trim_legend_title) {
            hlp_1$title <- gsub("\n.*",
               "",
               hlp_1$title)
         }
         hlp_1;
      })
      names(heatmap_legend_param) <- names(nmatlist);
   }
   if (any(c("legend_width", "border", "direction", "title_position") %in% names(heatmap_legend_param)) ||
         length(heatmap_legend_param) != length(nmatlist)) {
      if (any(c("legend_width", "border", "direction", "title_position") %in% names(heatmap_legend_param))) {
         heatmap_legend_param <- rep(list(heatmap_legend_param),
            length.out=length(nmatlist));
      } else {
         heatmap_legend_param <- rep(heatmap_legend_param,
            length.out=length(nmatlist));
      }
   }

   ## Handle top_annotation provided as a custom list
   top_annotation_list <- NULL;
   if (is.list(top_annotation)) {
      top_annotation_list <- top_annotation;
      if (!all(names(nmatlist) %in% names(top_annotation_list))) {
         if (length(top_annotation_list) != length(nmatlist)) {
            top_annotation_list <- rep(top_annotation_list,
               length.out=length(nmatlist));
            names(top_annotation_list) <- names(nmatlist);
         }
      }
      top_annotation_list <- top_annotation_list[names(nmatlist)];
   }
   top_axis <- rep(TRUE, length(nmatlist));
   if (length(panel_groups) > 0) {
      if ("left" %in% top_axis_side) {
         top_axis <- c("blahblah", head(panel_groups, -1)) != panel_groups;
         top_axis_side <- ifelse(top_axis, "left", "left");
      } else if ("right" %in% top_axis_side) {
         top_axis <- (panel_groups != tail(c(panel_groups, "blahblah"), -1))
         top_axis_side <- ifelse(top_axis, "right", "left");
      } else if ("both" %in% top_axis_side) {
         top_axis_left <- c("blahblah", head(panel_groups, -1)) != panel_groups;
         top_axis_right <- (panel_groups != tail(c(panel_groups, "blahblah"), -1));
         top_axis <- (top_axis_left | top_axis_right);
         top_axis_side <- ifelse(top_axis_right, "right",
            ifelse(top_axis_left, "left", "left"));
      } else if ("none" %in% top_axis_side) {
         top_axis <- rep(FALSE, length.out=length(panel_groups));
         top_axis_side <- rep("left", length.out=length(panel_groups));
      } else {
         top_axis <- rep(TRUE, length.out=length(panel_groups));
         top_axis_side <- rep("right", length.out=length(panel_groups));
      }
   } else {
      if (length(top_axis_side) == 0) {
         top_axis_side <- "left";
      }
      top_axis_side <- rep(top_axis_side,
         length.out=length(nmatlist));
   }

   #############################
   ## Iterate each heatmap
   if (verbose) {
      jamba::printDebug("nmatlist2heatmaps(): ",
         "Iterating each heatmap.");
   }

   signal_name_hash <- list();
   top_legend <- list();

   EH_l <- lapply(seq_along(nmatlist), function(i){
      nmat <- nmatlist[[i]][rows,,drop=FALSE];
      signal_name <- attr(nmat, "signal_name");
      target_name <- attr(nmat, "target_name");
      ## Check for duplicated signal_name
      if (signal_name %in% names(signal_name_hash)) {
         # append Excel-style column alphabetic lettering
         signal_name_hash[[signal_name]] <<- 1;
         signal_name <- paste0(signal_name, "_rep",
            jamba::colNum2excelName(i))
      }
      signal_name_hash[[signal_name]] <<- 1;
      s_name <- gsub("_at_", "\nat_", signal_name);
      if (length(nmat_names) > 0) {
         signal_name <- nmat_names[i];
      } else {
         # otherwise consider adding transform_label suffix
         if (nchar(transform_label[[i]]) > 0 &&
               !any(c(NA, "none") %in% transform_label[[i]])) {
            signal_name <- paste0(signal_name, "\n",
               transform_label[[i]])
         } else if (any(nchar(transform_label) > 0)) {
            signal_name <- paste0(signal_name, "\n ")
         }
      }

      color <- nmat_colors[[i]];
      # 0.0.76.900
      # - Check if color is empty, or non-function with NA, then use default.
      # - Also change default to `"Reds"`.
      if (length(color) == 0 || (!is.function(color) && any(is.na(color)))) {
         # color <- "aquamarine4";
         color <- "Reds";
      }
      ## Define ylim
      if (length(ylims) > 0) {
         ylim <- sort(ylims[[i]]);
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
            "signal_name:\n'",
            signal_name,
            "',\ntarget_name:",
            target_name,
            ",\ncolor:",
            color_txt,
            ",\nylim=(", jamba::cPaste(ylim_txt), ")",
            fgText=c(
               rep(list("darkorange",
                  "dodgerblue"),
                  length.out=6),
               NA),
            bgText=as.list(
               rep(list(NA),
                  length.out=6),
               color_txt)
         );
      }
      if (length(iceiling) > 0 && !is.na(iceiling)) {
         iceiling <- get_nmat_ceiling(imat,
            iceiling,
            verbose=verbose>1);
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

      ## if partition is a factor, call factor() which forces
      ## it to drop any missing factor levels
      if (is.factor(partition[rows])) {
         partition <- factor(partition[rows]);
      }
      use_colors <- k_colors[levels(partition)];
      if (length(use_colors) == 0) {
         use_colors <- k_colors[unique(partition)];
         if (length(use_colors) == 0) {
            n <- max(c(1, length(unique(partition))));
            use_colors <- colorjam::rainbowJam(n);
         }
      }
      if (length(use_colors) == 1 && length(names(use_colors)) == 0) {
         names(use_colors) <- "Overall";
      }

      ######################################
      # process top_annotation
      if (length(top_annotation_list) > 0) {
         top_annotation <- top_annotation_list[[i]];
      }
      if (length(top_annotation) > 0 && is.logical(top_annotation) && !top_annotation) {
         top_annotation <- NULL;
      } else {
         if (length(profile_linetype) == 0) {
            profile_linetype <- 1;
         }
         if (length(profile_linewidth) == 0) {
            profile_linewidth <- 1.5;
         }
         # Define top_legend
         use_linewidths <- rep(profile_linewidth, length.out=length(use_colors))
         use_linetypes <- rep(profile_linetype, length.out=length(use_colors))
         # This assignment looks bad, but top_legend was defined within
         # this function, which means <<- will assign the value to the
         # internal function variable and not to .GlobalEnv
         top_legend <<- ComplexHeatmap::Legend(title="Profiles",
            title_position="topleft",
            type="lines",
            at=names(use_colors),
            labels=names(use_colors),
            border=FALSE,
            nrow=length(use_colors),
            grid_width=grid::unit(14, "mm"),
            background="white",
            title_gp=grid::gpar(
               fontsize=legend_fontsize[1],
               fontface="bold"),
            labels_gp=grid::gpar(
               fontsize=legend_fontsize[2]),
            legend_gp=grid::gpar(
               col=use_colors,
               lty=use_linetypes,
               lwd=use_linewidths)
         )
         # Define top_annotation
         top_annotation <- ComplexHeatmap::HeatmapAnnotation(
            lines=EnrichedHeatmap::anno_enriched(
               gp=grid::gpar(col=use_colors,
                  lty=profile_linetype,
                  lwd=profile_linewidth),
               value=profile_value,
               ylim=ylim,
               axis=top_axis[i],
               axis_param=list(side=top_axis_side[i],
                  gp=axis_name_gp[[i]]),
               height=top_anno_height,
               show_error=show_error)
         )
      }

      ######################################
      # create coverage heatmap
      if (FALSE %in% row_order) {
         row_order <- NULL;
      }
      EH <- EnrichedHeatmap::EnrichedHeatmap(imat[rows,],
         row_split=partition[rows],
         pos_line=pos_line[[i]],
         use_raster=use_raster,
         raster_quality=raster_quality,
         raster_by_magick=raster_by_magick,
         col=colramp,
         border=border[[i]],
         top_annotation=top_annotation,
         show_heatmap_legend=show_heatmap_legend[[i]],
         heatmap_legend_param=heatmap_legend_param[[i]],
         axis_name_gp=axis_name_gp[[i]],
         axis_name=axis_name[[i]],
         axis_name_rot=axis_name_rot,
         name=signal_name,
         column_title=signal_name,
         row_order=row_order[rows],
         row_title_rot=row_title_rot,
         column_title_gp=column_title_gp[[i]],
         ...);
      EH;
   });

   ######################################
   ## Create Summary as Custom Legend
   if (TRUE) {
      use_main_heatmap <- main_heatmap;
      if (all(seq_along(nmatlist) %in% main_heatmap)) {
         use_main_heatmap <- " (All)";
      } else {
         use_main_heatmap <- jamba::cPasteSU(main_heatmap, sep=", ")
      }
      caption_list <- list(
         Summary=jamba::rmNA(c(
            ifelse(length(nmatlist) > 1,
               paste0("Heatmaps: ", length(nmatlist)),
               NA),
            paste0("Rows: ",
               jamba::formatInt(nrow(nmatlist[[1]]))),
            ifelse(length(nmatlist) > 1,
               paste0("Main Heatmap",
                  ifelse(length(main_heatmap) > 1, "s:\n    ", ": "),
                  use_main_heatmap),
               NA)
         ))
      )

      # Optional row sorting
      if (length(byCols) > 0) {
         caption_list$Summary <- c(
            caption_list$Summary,
            c(
               paste0("Sorted By:\n    ",
                  jamba::cPaste(byCols, sep=",\n    "))
            )
         )
      }

      # describe partitioning when used
      if (length(partition) > 0 && length(unique(partition)) > 1) {
         caption_list$Partitions <- c(
            c(
               paste0("Partitions: ",
                  length(unique(partition)))
            )
         )
      }
      if (any(k_clusters > 1)) {
         caption_list$Partitions <- c(
            caption_list$Partitions,
            c(
               paste0("k-Method: ", k_method),
               paste0("k-Heatmap",
                  ifelse(length(k_heatmap) > 1, "s:\n    ", ": "),
                  jamba::cPasteSU(k_heatmap))
            )
         )
      }
      # describe recentering and restranding when used
      if (length(recenter_heatmap) > 0) {
         caption_list$Adjustments <- c(
            caption_list$Adjustments,
            paste0("Re-Center Heatmap",
               ifelse(length(recenter_heatmap) > 1, "s:\n    ", ": "),
               jamba::cPasteSU(recenter_heatmap, sep=", "))
         )
      }
      if (length(restrand_heatmap) > 0) {
         caption_list$Adjustments <- c(
            caption_list$Adjustments,
            paste0("Re-Strand Heatmap",
               ifelse(length(restrand_heatmap) > 1, "s:\n    ", ": "),
               jamba::cPasteSU(restrand_heatmap, sep=", "))
         )
      }

      # make convenient text summary
      caption <- paste0(collapse="\n",
         paste0(names(caption_list), ":\n  "),
         jamba::cPaste(caption_list, sep="\n  "));
      # Convert to list of Legend objects
      caption_legends <- lapply(names(caption_list), function(lname){
         ilabel1 <- caption_list[[lname]];
         ilabel <- jamba::cPaste(sep="\n",
            caption_list[[lname]]);
         caption_grob <- grid::textGrob(
            label=jamba::cPaste(sep="\n",
               caption_list[[lname]]),
            gp=grid::gpar(
               fontsize=legend_fontsize[2],
               lineheight=1),
            hjust=0, vjust=0)
         text_legend2 <- ComplexHeatmap::Legend(
            title_gp=grid::gpar(
               fontsize=legend_fontsize[1],
               fontface="bold"),
            title=paste0(lname),
            grob=caption_grob,
            title_position="topleft",
            type="grid")
      })
      # Convert to LegendList
      caption_legendlist <- ComplexHeatmap::packLegend(
         list=c(caption_legends),
         direction="vertical")
   }

   ######################################
   ## Layout heatmap panels
   ##
   HM_drawn <- NULL;
   if (length(hm_nrow) > 0 &&
         hm_nrow > 1 &&
         length(nmatlist) > 1) {
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

   # 0.0.76.900 - store important function parameters in returned object
   anno_df_colnames <- NULL;
   if (length(anno_df) > 0) {
      anno_df_colnames <- colnames(anno_df);
   }
   #
   fn_params <- list();
   fn_params$rows <- rows;
   fn_params$anno_df_colnames <- anno_df_colnames;
   fn_params$partition <- partition;
   fn_params$k_clusters <- k_clusters;
   fn_params$k_colors <- k_colors;
   fn_params$nmat_colors <- nmat_colors;
   fn_params$transform <- transform;
   fn_params$panel_groups <- panel_groups;
   fn_params$ylims <- ylims;
   fn_params$signal_ceiling <- signal_ceiling;
   fn_params$ht_gap <- ht_gap;
   fn_params$top_annotation_list <- top_annotation_list;
   if (length(recenter_heatmap) > 0) {
      fn_params$recenter_heatmap <- recenter_heatmap;
      fn_params$recenter_range <- recenter_range;
      fn_params$recenter_invert <- recenter_invert;
   }
   if (length(restrand_heatmap) > 0) {
      fn_params$restrand_heatmap <- restrand_heatmap;
      fn_params$restrand_range <- restrand_range;
      fn_params$restrand_buffer <- restrand_buffer;
      fn_params$restrand_invert <- restrand_invert;
   }

   # 0.0.76.900 - include visual caption summary
   # - number of rows, k-means method
   hm_caption <- caption;
   if (FALSE) {
      hm_caption <- paste0(
         jamba::formatInt(length(nmatlist)), " heatmap",
         ifelse(length(nmatlist) > 1, "s", ""), ", ",
         jamba::formatInt(length(rows)), " rows")
      if (length(partition) > 0) {
         np <- length(unique(partition));
         hm_caption <- paste0(hm_caption,
            "\n", np,
            " partition", ifelse(np > 1, "s", ""));
      }
      if (length(byCols) > 0) {
         hm_caption <- paste0(hm_caption,
            "\nsorted by:\n   ",
            jamba::cPaste(byCols, sep="\n   "))
      }
      if (any(k_clusters > 1)) {
         hm_caption <- paste0(hm_caption,
            "\nk-means '", k_method, "' ",
            "\n   using heatmap", ifelse(length(k_heatmap) > 1, "s", ""), " ",
            jamba::cPasteS(k_heatmap))
      }
      hm_caption <- paste0(hm_caption,
         "\nmain heatmap", ifelse(length(main_heatmap) > 1, "s", ""),
         " ", jamba::cPasteS(main_heatmap))
      # draw_function()
      draw_caption <- function
      (text=hm_caption,
       x=grid::unit(1, "npc"),
       y=grid::unit(0, "npc"),
       fontsize=10,
       font="Arial",
       just=c("right", "bottom"),
       color="midnightblue",
       background_color="white",
       ...)
      {
         drawn <- suppressWarnings({
            ComplexHeatmap::grid.textbox(
               text=text,
               x=x,
               y=y,
               just=just,
               gp=grid::gpar(
                  col=color,
                  fontsize=fontsize,
                  font=font),
               background_gp=grid::gpar(
                  col="transparent",
                  fill=background_color),
               ...)
         })
         return(invisible(drawn))
      }
      if (TRUE %in% do_plot && TRUE %in% do_caption) {
         draw_caption(fontsize=legend_fontsize[2]);
      }
   }

   ret_list <- list(
      AHM=AHM,
      PHM=PHM,
      EH_l=EH_l,
      MHM=MHM,
      draw=list(
         HM_temp=HM_temp,
         ht_gap=ht_gap,
         main_heatmap=main_heatmap_temp),
      caption_legendlist=caption_legendlist
   );
   # return HM_drawn with the heatmap as drawn
   if (length(HM_drawn) > 0) {
      ret_list$HM_drawn <- HM_drawn;
   }
   ret_list$fn_params <- fn_params;
   ret_list$hm_caption <- caption;
   ret_list$adjust_df <- adjust_df;
   # ret_list$draw_caption <- draw_caption;
   invisible(ret_list);
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
#' @returns `function` or `NULL` when no matching function is
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
      if (any(c("", NA, "none") %in% it)) {
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
#' It takes a `normalizedMatrix` or `numeric` matrix object, and
#' a ceiling value `iceiling` and determines an appropriate numeric
#' ceiling with the following rules:
#'
#' * if `iceiling` is `NULL`, it returns the highest absolute value in `imat`
#' * if `iceiling > 0` and `iceiling <= 1`, it calculates
#' `quantile(abs(imat), probs=iceiline)`, excluding values of zero
#' * otherwise `iceiling` is used as a fixed numerical ceiling
#'
#' In all cases, `iceiling` is rounded to 3 digits with `round(iceiling, digits=3)`
#'
#' Also in all cases, `na.rm=TRUE` is used, to prevent returning `NA`.
#'
#'
#' @family jam coverage heatmap functions
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
         jamba::printDebug("get_nmat_ceiling(): ",
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
      jamba::printDebug("get_nmat_ceiling(): ",
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
