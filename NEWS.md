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
