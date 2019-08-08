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
