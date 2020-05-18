# platjam version 0.0.26.900

## changes to existing functions

* `nmatlist2heatmaps()` now treats annotation columns in `anno_df`
as bi-directional when the max value is <= 50, which helps in cases
that have values `c(-1, 0, 1)` and sometimes only `c(0, 1)`.
The end result is more numeric columns will have the same color
gradient (blue-white-red).
* `nmatlist2heatmaps()` displays legend with distinct color steps
for `anno_df` columns with 10 or fewer distinct values, which is
mostly helpful with a small number of integer values.

# platjam version 0.0.25.900

## bug fixes

* `applyXlsxConditionalFormatByColumn()` fixed bug in verbose
output that referred to nonexistent object. Also removed `verbose=TRUE`
from calling function `save_salmon_qc_xlsx()`.

# platjam version 0.0.24.900

## new functions

* `save_salmon_qc_xlsx()` saves output from `get_salmon_meta()`
in a convenient Excel `.xlsx` format, with four worksheets showing
different aspects of data quality. It returns a `list` of
`data.frame` summaries invisibly, for further review.
* `applyXlsxConditionalFormatByColumn()` applies conditional
formatting using values in each column, rather than using
a fixed range for all columns. This function will become
part of `jamba::writeOpenxlsx()` and
`jamba::applyXlsxConditionalFormat()` in the near future.
* `set_xlsx_colwidths()` and `set_xlsx_rowheights()` define
fixed column width and row height, respectively, for `.xlsx` files.

# platjam version 0.0.23.900

## new functions

* `color_complement()` takes input color, rotates the color hue 180
degrees to make a complement color. It is unique in that it uses
warped hue (see `colorjam::h2hw()`) which therefore uses red-yellow-blue
color wheel instead of red-green-blue, which is much better for
color operations like these. This function will almost certainly
be moved into the `colorjam` package. The examples show how to
make a divergent color gradient from a single color.

## changes to existing functions

* `nmatlist2heatmaps()` new argument `nmat_names` used to supply
custom names for each heatmap, otherwise the attribute
`"signal_name"` values are used.
* `nmatlist2heatmaps()` fixed error caused by `use_raster=TRUE`
and any operation that expects to create a grid graphical object,
such as `multienrichjam::grid_with_title()`, or multi-row output
specified with `hm_nrow=2`.
* `nmatlist2heatmaps()` new arguments `heatmap_legend_param` and
`annotation_legend_param` allow customizing the legends in detail,
passed to `ComplexHeatmap::Heatmap()`. However, default settings
show color gradients horizontally, which may save some screen space.
* `nmatlist2heatmaps()` new arguments `title` and `caption` invoke
a prototype method that calls `multienrichjam::grid_with_title()`
to position an overall title and caption above and below the
heatmaps, respectively. When supplied, these sections are displayed
and the heatmaps are resized accordingly.
* `nmatlist2heatmaps()` now properly handles `k_subset` to plot
a subset of clusters provided by `partition` or after `k_clusters`
is used to create clusters. Mostly useful for follow-up plots
to drill down into one or more specific clusters.
* `nmatlist2heatmaps()` now properly handles `row_order` for edge
cases like `k_subset`, or when `row_order=NULL`; also `row_order=TRUE`
will enable the default `EnrichedHeatmap::enriched_score()` for
the `main_heatmap` using only the rows being displayed.
* `nmatlist2heatmaps()` new argument `border` whether to draw a border
around the heatmap. It can be supplied as a vector, to control the
border around each heatmap if needed.

# platjam version 0.0.22.900

## changes to existing functions

* `nmatlist2heatmaps()` accepts input `nmat_colors` with
color functions, for example sending output from `circlize::colorRamp2()`.
* `nmatlist2heatmaps()` when given a single color per matrix,
when data is determined to have positive and negative values,
it generates a divergent color ramp using the supplied color,
and complementary color. To define a specific set of colors,
supply a three color vector, the middle color will be fixed
at zero.

## bug fixes

* `nmatlist2heatmaps()` no longer requires `ylims` and uses `ylim=NULL`
properly.

# platjam version 0.0.21.900

## changes to existing functions

* `parse_ucsc_gokey()` was updated to recognize bigBed/bb files
can not become overlay tracks.

# platjam version 0.0.20.900

## changes to existing functions

* `get_track_defaults()` was modified to include `shortLabel,longLabel`
in the template for multiWig tracks. Apparently a label is displayed
and it uses the internal hub number as a prefix by default.

# platjam version 0.0.19.900

## new functions

* `get_salmon_meta()` is intended to load Salmon transcript
quantification metadata from several output JSON files, returning
results in a `data.frame` with one row per Salmon run. It combines
various measured statistics from `"meta_info.json"`, `"cmd_info.json"`,
and `"lib_format_counts.json"` so that multiple Salmon runs can be
directly compared in a table format.
* `get_salmon_root()` is a wrapper to `rprojroot::find_root()`, intended
to find the root directory of a Salmon run given any Salmon output
file. It searches for the root directory that contains the file
`"cmd_info.json"`.

# platjam version 0.0.18.900

## Changes to existing functions

* `get_numeric_transform()` was updated to clean up the help
docs, to accept `"linear"` as a valid name (which does not change
the input data), and to parse the input `transform` values in
order to assign reasonable names where possible.

# platjam version 0.0.17.900

## Bug fixes

* Fixed various bugs in `nmatlist2heatmaps()`. No doubt there will
be more, covering the host of assumptions made for the panel_groups,
signal_ceiling, and ylims arguments.

# platjam version 0.0.16.900

## Changes to `nmatlist2heatmaps()`

* Cool new feature: `nmatlist2heatmaps()` argument `panel_group`,
when it is defined, it calculates a shared `ylims` for each group
of panels, and a shared `signal_ceiling` for each group of panels.
The signal ceiling uses the supplied `signal_ceiling` for each panel
(if supplied) then takes the maximum value for each group of panels.
So if `signal_ceiling=0.8` it is treated as quantile 0.8, the 0.8
quantile is calculate for all non-zero values for each panel, then
each group of panels uses the max for that group.
* `nmatlist2heatmaps()` the `ht_gap` gap between heatmaps is not
applied between annotation heatmaps, instead `"1mm"` gap is used
between annotation heatmaps. The purpose of `ht_gap` is to have
enough gap between heatmaps so the y-axis profile plot labels
can be displayed, but there are no y-axis labels for annotations.
* The default `lens` for `anno_df` numeric columns is now `lens=1`,
reverting previous change `lens=2` which was causing gradients to
be darker than expected.
* The verbose output is now more concise, and TRUE by default. For
even more detail use `verbose=2`.

## New functions

* `get_nmat_ceiling()` is a simple wrapper function called by
`nmatlist2heatmaps()` and not intended to be called directly.
It takes a numeric matrix, and ceiling, and returns what I determines
to be the appropriate ceiling to use.

# platjam version 0.0.15.900

## Changes to existing functions

* `nmatlist2heatmaps()` modified `signal_name` to clean up leading/trailing
whitespace and return `\n`, and repeated return `\n` characters.
* `nmatlist2heatmaps()` new argument `ht_gap` used to define the gap between
heatmap panels. This argument is intended to help reduce the overlap of
the top y-axis profile labels in adjacent profile panels.
* `nmatlist2heatmaps()` new argument `axis_name`, as either a vector, or
a list of vectors, used for x-axis labels across each panel.
* `nmatlist2heatmaps()` new argument `column_title_gp`, as either a
`grid::gpar` object, or list of gpar objects.
* `nmatlist2heatmaps()` new argument `profile_value`, passed to
`EnrichedHeatmap::anno_enriched()` to define the type of numeric summary
to use for the profile plots.
* `nmatlist2heatmaps()` new argument `panel_groups`, character vector of
groups used to apply consistent y-axis ylims values per group. By default
if `nmat_colors` is not supplied, then `colorjam::group2colors(panel_groups)`
is used to define group colors.

# platjam version 0.0.14.900

## Bug fixes

* `nmatlist2heatmaps()` fixed bug in handling argument `signal_ceiling`,
which when the value is between 0 and 1 is treated as a quantile. The
method now removes zero before the `quantile()` calculation.
* `nmatlist2heatmaps()` fixed bug that wrongly expanded argument
`axis_names_gp` as a list without realizing class `"gpar"` is also
a list. Now `"gpar"` is wrapped in `list()` before expanding to the
number of coverage matrices.
* `nmatlist2heatmaps()` fixed bug where argument `axis_name_rot` was
not defined in the function itself.

# platjam version 0.0.13.900

## Enhancements

* `nmatlist2heatmaps()` altered the logic used to define color ranges
in numeric columns of `anno_df`, opting to use quantile ranges
`c(0.005, 0.995)` to help trim extreme values. Also changed `"purple"`
to `"Purples"` to use the `RColorBrewer` purple color gradient.

# platjam version 0.0.12.900

## bug fixes

* `get_numeric_transform()` bug fixed wrong argument name used inside
the function. Should also resolve issues with `nmatlist2heatmaps()`.

# platjam version 0.0.11.900

## new functions

* `import_nanostring_rcc()` reads Nanostring RCC files, and produces
a `SummarizedExperiment` object suitable for use in omics-style
analyses.

# platjam version 0.0.10.900

## enhancements

* `nmatlist2heatmaps()` new argument `signal_ceiling` defines
a numerical ceiling for each input matrix, used to apply the
color ramp to a defined numeric range. When `signal_ceiling`
is between 0 and 1, it is assumed to be a quantile, so the
appropriate quantile value is used as the ceiling. These
values can be seen when `verbose=TRUE`.
* `nmatlist2heatmaps()` argument `transform` now accepts text
string for transformation functions, such as `"log2signed"`,
`"sqrt"`, `"cubert"`, etc. Any other string will try to
match a function name in the R search path.

## new functions

* `get_numeric_transform()` is used to apply name-based
transformations, called from `nmatlist2heatmaps()`.

## UCSC genome browser functions

Some new functions were added to help create UCSC genome
browser track hubs, which are fairly tedious to build manually.

* `get_track_defaults()` returns default values for tracks,
including some track snippet templates for composite and
overlay track types.
* `parse_ucsc_gokey()` is a parser for a specific track type,
based upon prominent scientist Dr. Gokey, who defines a straightforward
format for a series of grouped tracks.
* `make_ucsc_trackname()` creates a valid UCSC track name, removing
invalid characters.
* `assign_track_defaults()` assigns track values to an environment
which may already contain default values. New values override previous
values, after which the environment can be used with the `glue`
package to update a track snippet template for each track in a set.


# platjam version 0.0.9.900

## enhancements

* `coverage_matrix2nmat()` now accepts a vector of files as
input, and will process them and return a list of normalizedMatrix
objects, ready to go directly into `nmatlist2heatmaps()`.

# platjam version 0.0.8.900

## bug fixes

* `nmatlist2heatmaps()` fixed bug in placement of argument
show_error, now properly placed inside `EnrichedHeatmap::anno_enriched()`.

# platjam version 0.0.7.900

## changes

* `nmatlist2heatmaps()` new argument `show_error` is
passed to `EnrichedHeatmap::anno_enriched()` which enables
shaded error bars using standard error around the profile
at the top of each heatmap.

## new functions

* `frequency_matrix2nmat()` converts a frequency matrix
whose colnames represent numeric frequency values, and
converts to class `"normalizedMatrix"` for use in
EnrichedHeatmap. It optionally scales each row using
`rowNormScale()`.
* `rowNormScale()` applies `jamba::normScale()` to each
row of a matrix.

# platjam version 0.0.6.900

## changes

* `nmatlist2heatmaps()` argument `k_subset` now
works when `partition` is supplied and k-means
clustering is not performed inside this function.
Therefore, one can supply their own clusters, and
choose a subset for subsequent drilldown. The
colors are assigned by total cluster, so the subset
should retain consistent colors per cluster.
* `nmatlist2heatmaps()` argument `k_method` allows
defining different k-means clustering methods, including
`c("pearson","correlation")` if the `"amap"`
package is installed.

# platjam version 0.0.5.900

## changes

* Version bump while debugging logistics with 
`nmatlist2heatmaps()` argument `transform`.

# platjam version 0.0.4.900

## changes

* `nmatlist2heatmaps()` argument `transform` is now applied
to each heatmap, with `transform` values recycled to the
`length(nmatlist)`, to allow control over mathematical
transformation of each numeric matric.

# platjam version 0.0.3.900

## changes

* `nmatlist2heatmaps()` argument anno_df is a data.frame with
annotations that can be used to order rows in the data.
* `nmatlist2heatmaps()` updated its method to handle proper
order of user-supplied partition information.
* `nmatlist2heatmaps()` new arguments to handle row labeling,
still testing the combination of partitions and row labels
which appear misaligned.

# platjam version 0.0.2.900

## changes

* The package itself allows installation on R-3.0.0 and higher,
resolving an issue with installing on R-3.5.0 when the
package said it required R-3.5.1 (and did not).

# platjam version 0.0.2.900

## new function

* `nmatlist2heatmaps()` creates multiple coverage heatmaps
as the result of importing data with `coverage_matrix2nmat()`.

# platjam version 0.0.1.900

Initial release.

## new functions

* `coverage_matrix2nmat()` imports genome coverage matrix
file and converts to class `normalizedMatrix` for use in
the amazing `"EnrichedHeatmap"` package.
