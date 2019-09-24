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
