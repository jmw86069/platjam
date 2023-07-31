# TODO 31jul2023

* `design2colors()`

   * DONE. accept `DataFrame` input, as from `SummarizedExperiment::colData()`.
   * DONE. use workaround for class/group cardinality different than 1-to-N
   * DONE. when lightness is defined and group/lightness cardinality is not
   1-to-N, assign "grey45" to lightness_colnames in color_sub.
   * DONE. handle using `"-"` prefix for colnames to reverse the sort order
   * DONE. fix use case with no values assigned from `group_colnames`
   * fix use case with only one input colname, no rownames are printed
   * consider sorting column assignment based upon cardinality with
   `group_colnames`, then number of unique entries, in order to impose
   some reproducibility even when input column order is changed.
   * accept color assignment for `class_colnames` values, which would
   then impose constraints on group color hues within reasonable range
   of those colors. It would effective split the color hue into a slightly
   wider range of hues, centered on the class color hue.

# TODO 19jul2023

* DONE: add importer for metabolomics data

   * `import_metabolomics_niehs()` - import LC-MS or LC-MS/MS data
   processed by the NIEHS Metabolomics Core Facility, with defined file
   formats.

# TODO 18apr2023

* `nmatlist2heatmaps()`

   * add method to display group labels above heatmaps:
   
      * `signal` - (H3K27ac, H3K4me3, etc) across contiguous `panel_groups`
      * `type` - (signal, difference) across contiguous signal types
      within `panel_groups`
      * `label` - trimmed label which no longer requires the `signal` and
      `type` encoded into the name.
      * use `heatmap_column_group_labels()` style
      
   * add method to adjust padding (whitespace) around the overall heatmap
   layout; between heatmaps and color legend; between annotation stripes.
   * make the color assignment more similar to `jamses::heatmap_se()`,
   in the form of color `color_list` named by column.

# TODO 22feb2023

* `design2colors()`

   * Error is thrown when there are empty values in `group_colnames` fields.
   * DONE: Allow using color ramp name instead of color for column assignment.
   * DONE: Allow defining a `color_max`, maximum value for `numeric` columns that
   will use a color gradient.
   
      * This option will help whenever a few columns should share the
      same threshold. It is not auto-detected, but can still be useful.
      The analyst will need to provide the `color_max` value upfront,
      perhaps in the form of a named vector similar to `color_sub`,
      named by `colnames(x)`.
   
   * Allow `"rownames"` as valid input to `group_colnames`,
   essentially assigning one value per row.
   * Make the function work when `group_colnames` is `NULL`.
   Assign colors to rownames?
   * consider removing `class_group_color` and `class_group_lightness_color`
   from output, or renaming each to use the combination of actual colnames
   assembled to form those colors.

# TODO 10jan2023

* `get_salmon_meta()` - please fix the file ordering as described below.

# TODO 20sep2022

* `get_salmon_meta()` should return data in the same order the filenames
were provided, in fact it's weird that it would not already do that.

# TODO 13sep2022

* (COMPLETE) Add import function for Mascot proteomics type data format.

# TODO 13may2022

* `design2colors()`

   * make it work even without `group_colnames`, which is effectively
   just colorizing a `data.frame` by each column.
   * Make `numeric` columns apply a gradient as if a categorical
   color were supplied. Basic idea is that instead of assigning
   categorical colors to each numeric value, for numeric columns
   assign the categorical color to the column name. Then it should
   proceed as it does now, when the user supplies a `color_sub`
   for that column name.

# TODO 05apr2022

## Salmon import to SummarizedExperiment?

* goal would be to mimic `import_proteomics_PD()` by importing
`quant.sf` files, incorporate `curation.txt` for sample annotations.

   * requires `GTF` and/or `tx2gene` file for transcript-to-gene.
   * Produce `list` with `TxSE` and `GeneSE` objects.
   * `assays()` will include `counts` and `abundance`.

## Salmon metadata

COMPLETE: Salmon produces a useful file `"flenDist.txt"` that includes the
distribution of fragment lengths observed, from 1 to 1000 length.
Goal is to parse this file to determine the weighted mean
fragment length, likely extending `get_salmon_meta()`.


# TODO 21mar2022

## UCSC track hubs

* `parse_ucsc_gokey()` updates:

   * currently the compositeTrack output requires editing to become
   visible by default. It appears now to export a mix of
   superTrack/compositeTrack. The apparent changes required:
   
      1. Remove "compositeTrack on" from the 2nd level of compositeTrack view.
      2. Remove "superTrack on show" from the 1st level of compositeTrack view.
      3. Add "type bigWig" to the 1st level of compositeTrack view (unlikely 
      to be necessary however).
      4. Add to 3rd level compositeTrack view, for each track:
      `"subGroups            view=COV"`
      Unclear if this line is required for visibility.

   * Option `autoScale group` is useful for sets of tracks, and it should
   be easy to configure upfront.
   
      * Track groups should be possible to be created independent of
      track order. For example ATAC_cov and ATAC_NDR could appear for
      each sample in order, with ATAC_cov autoScale to its own group,
      and ATAC_NDR autoScale to its own group.
   
   * Not all track arguments are being honored in the final track hub entry.
   For example autoScale=group is not propagating through the workflow.



# TODO for platjam

## `nmatlist2heatmaps()`

* when `transform` is used, alter the y-axis numeric labels,
and color legend numeric labels accordingly.
* Optionally draw y-axis horizontal abline at y=0
when data contains positive/negative values.


## UCSC track hub

* Add example for `parse_ucsc_gokey()` showing three kinds of tracks
as input

   1. bigWig with pos/neg or F/R or some track grouping -> multiWig overlay
   2. bigWig without grouping -> composite tracks
   3. bigBed -> composite Tracks

* COMPLETE:new arguments for bigBed tracks

   * scoreFilter=0, scoreFilterLimits="0:1000"
   * minGrayLevel=4
   * spectrum="off", scoreMin=0, scoreMax=1000


# TODO for platjam 28aug2020


## Specific usability

* For `panel_groups`:

   * COMPLETE: when all heatmaps in a set of panels share the
   same color gradient, the color legend should be displayed once per group,
   and should be labeled using the `panel_groups` value.
   * make the `ht_gap` between adjacent heatmaps in the same group smaller
   * COMPLETE: hide y-axis labels on the profile plot for all but the last panel in
   an adjacent set. For example using
   `EnrichedHeatmap::anno_enriched(... yaxis=FALSE)`


## Visual enhancements

* `anno_df` should probably use `Heatmap()` and not `HeatmapAnnotation()`
so that it can enable `use_raster=TRUE`.
* Color legend for the profile plot should use lines instead of boxes
to indicate colors (if possible).


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
