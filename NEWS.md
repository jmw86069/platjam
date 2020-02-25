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
