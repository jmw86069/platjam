# TODO for platjam 28aug2020

## k-means and partitioning, correct the order

Currently k-means clustering and partitioning works together,
but the k-means is performed on everything, then clusters are
divided by partition.

The effect when partitioning by "up" and "down" is that some
k-means clusters are dominated by the "up" and "down" such
that cluster 1 has 95% members from "up", very few from down.
For k=4, it ends up making 8 clusters, which are mostly just
4 proper clusters, and 4 junk clusters (the "up" cluster with few
"down" members).

Would be better to k-means within each partition, to make proper
subclusters for each partition.

# TODO for platjam 22may2020

## Dynamic ylims (DONE)

* The detected ylim range should include the additional
range when show_error=TRUE. This change involves calculating
the error bars during the ylim calculations, not too
computationally intensive, but it must replicate the
method used by EnrichedHeatmap.

## Custom color substition (colorSub) (DONE)

* `nmatlist2heatmaps()` should allow optional `colorSub`
argument as a named vector of colors. Any categorical column
all of whose values match `names(colorSub)` will use colors
in `colorSub` instead of generating its own new categorical
colors. This option will help allow colors to be consistent.

## Allow both partition and kmeans clustering (DONE)

* Currently kmeans clustering overrides partition, but we
want the option to cluster and partition.

## Genome coverage matrix methods (DONE)

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
