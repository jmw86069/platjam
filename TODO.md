# TODO for platjam 28aug2020


## Big picture - make it easy to use

Make coverage heatmaps easier for other scientists.

* tab-delimited table input with the following columns:

   * coverage matrix file path - full path to coverage file
   * name - name associated to the coverage matrix file, used below
   * display - logical indicating whether to display this coverage
   matrix; TRUE, 1 means "yes"; FALSE, 0 means "no". When "no", the
   coverage may still be usable in difference matrix calculations.
   * label - need some way to specify line breaks, without using line breaks.
   * xaxis_names - three comma-delimited values, or `NA`, for example:
   `-250kb,center,+250kb`
   * group - used for `panel_groups`
   * signal_ceiling
   * ylim
   * color_gradient - string that matches jam_linear, jam_divergent,
   viridis, RColorBrewer, or a comma-delimited series of hex colors
   * control_name - optional name of a control coverage file to be used
   to create a difference matrix

* tab-delimited anno_df file whose first column values should match
rownames in each coverage matrix file.

Make jam color gradients available:

   * `jam_linear` - color blind friendly linear gradients
   from white to color
   * `jam_divergent` - color blind friendly divergent gradients
   with black as the central color
   * The gradients above are designed so that the divergent
   start or end color matches the color name of linear gradients.
   The intent is to use linear gradients to represent straight
   coverage, and divergent gradients for coverage differences.

Make it possible to define coverage differences

* Probably use some format in the tab-delimited file above.


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
