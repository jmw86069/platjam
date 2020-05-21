# TODO for platjam

## Genome coverage matrix methods

* `nmatlist2heatmaps()` argument for `data.frame` of annotations,
used to display alongside heatmaps, and/or used to sort the heatmap
rows. For example one column could contain `"log2foldchange"`,
or a categorical variable.
* add a new file importer for deepTools coverage matrix files.
Allow for multiple samples to be present in the same file, then
create one `normalizedMatrix` for each sample with consistent
rownames and colnames.

## Genome coverage profile methods

* Similar to `nmatlist2heatmaps()` except that it focuses on
just the profile plots, including optional error bars and
statistical testing by position.

## Custom x-axis range

* Ability to specify a custom subset range for the x-axis.
Note that matrices in `nmatlist` are not required to share
the same x-axis range, so this function would probably
need to be applied to individual `normalizedMatrix` entries.
