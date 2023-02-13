
# Ideas for coverage heatmap Rmarkdown

* user-defined output file
* option to use random subset of N regions
* subtraction of two signals (B-A)
* tabular input `coverage_file`

   * file (or contrast)
   * code (preferably a letter)
   * panel_group - `character` group name used to associate multiple
   coverage lanes with the same settings: `color_ramp`, `ylims`,
   `signal_ceiling`. By default, color legend is shown for only the
   first member of each `panel_group`.
   * color_ramp - comma-delimited colors, or name of color ramp recognized
   by `jamba::getColorRamp()`, which can be tested with
   `jamba::showColors(jamba::getColorRamp("RdBu_r"))`
   * transform - `character` name of transformation function:
   `"log2signed"`, `"sqrt"`, `"cubert"`, `"qrt"`.
   * signal_ceiling - `numeric` color maximum for the color ramp
   * axis_name - `character` comma-delimited labels for the x-axis, usually
   three labels. Default uses the data encoded in the coverage file range.
   * axis_name_fontsize - `numeric` fontsize in points
   * axis_name_rot - `numeric` rotation angle
   * column_fontsize - `numeric` fontsize in points
   * lens - `numeric` adjustment to color ramp sensitivity, default `lens=-2`
   * pos_line - 1,0 whether to draw line at the middle position
   * ht_gap - `numeric` gap after this heatmap in "mm"
   * profile_value - `"mean", "sum", "abs_mean", "abs_sum"`
   * show_heatmap_legend - 0,1 whether to display color legend for this panel.
   Default shows color legend for the first panel in each `panel_group` when
   defined.
   # * ylims
   # * border
   # * legend_width
   # * trim_legend_title
   # * heatmap_legend_param

   * yaxis range
   * xaxis labels
   * upstream_length,downstream_length - to override data ranges
   * top_axis_side

* other arguments

   * title - `character` overall title
   * figure_height,figure_width - `numeric` height, width in inches.
   Default will use `figure_height=8` cm, then width:

      * 0.8 inch per heatmap panel
      * 1.5 inch total color legends
      * 0.2 inch per value in `anno_colnames`
   
   * top_anno_height - `numeric` height in `cm`
   * legend_max_labels - `integer` maximum color legend labels, in case
   there are 40+ heatmaps
   * hm_nrow - `integer` number of heatmap layout rows
   * seed - `numeric` seed used to define "reproducible randomness"
   * use_raster - 0,1 whether to convert output to rasterized image,
   almost always necessary and beneficial
   * raster_quality - `integer` factor used to combine adjacent numeric
   matrix values prior to converting to image output. Useful when
   the R `"magick"` package is not available.
   * raster_by_magick - 0,1 whether to use the R `"magick"` package if
   available.
   * verbose - 0,1 whether to enable verbose output during processing,
   sometimes useful for troubleshooting.

* row partitioning arguments

   * k_clusters - `integer` number of k-means clusters, which therefore
   enables k-means clustering.
   * k_subset - `integer` number or comma-delimited numbers of k-means
   clusters to display after calling k-means clustering.
   This option is used as a follow-up step, for example running k-means once
   then re-running to zoom into one or more specific clusters.
   * k_colors - `character` comma-delimited colors to assign to k-means
   clusters, in the order they appear.
   * k_width - `numeric` width of k-means clustering annotation stripe,
   in "mm".
   * k_method - `character` k-means method to use, typically: `"euclidean"`,
   `"correlation"`, but any method recognized by `amap::hcluster()`.
   * k_heatmap - `integer` or comma-delimited integers of the heatmap(s)
   to use during k-means clustering. Typically one heatmap is used for k-means
   clustering, then the partitioning is applied across all heatmaps.
   * partition - `character` of colname(s) in the optional `annofile` to use
   for partitioning the heatmap rows.


* optional annotation file arguments

   * annofile - `filename` full path to tab-delimited annotation data,
   whose first column should contain rownames that match rownames in
   the coverage matrix files, with additional columns that can be displayed
   alongside the coverage heatmap panels.
   * anno_colnames - `character` colnames in `annofile` to be included
   in the left annotation beside coverage heatmap panels. When no
   `anno_colnames` are define, none are displayed... to prevent displaying
   150 annotations, leaving no room to display coverage heatmaps.

* optional color file

   * name:color association
