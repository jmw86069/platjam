# TODO 18apr2025

* `import_salmon_quant()`

   * Allow importing data without GTF file, for example "tx" only.
   * Add `SummarizedExperiment::SummarizedExperiment()` prefix. Ugh.
   * Debug `curation_txt` argument when provided as filename. Ugh.

# TODO 24mar2025

Grand ideas for `rmd_tab_iterator()`:

* Create battery of functions to use with  for
Rmarkdown Omics analysis workflows:

   * `maplot_tab_fn()`
   * `heatmap_tab_fn()`
   * `volcanoplot_tab_fn()`
   * `covheatmap_tab_fn()`
   * Each could provide optional list to populate with data for inspection.

* Consider method to use `knitr::knit_child()` to allow custom output
figure sizes.
See: https://github.com/quarto-dev/quarto-cli/discussions/11042
* Try to silence multi-level tabs when all child tabs are also empty.
(This seems difficult*.)

* Improve `design2colors()`

   * It should be able to "figure out" a reasonably good default from
   inspecting the data.
   * It should more intuitively honor the given settings.
   * Add common defaults: Male/Female (M/F), Treatment, Time, Age,
   Control/Vehicle/Untreated, etc.

# TODO 21feb2025

* Consider `heatmaps2nmatlist()` as an inverse to `nmatlist2heatmaps()`
* Review ways to cache/reuse the recentering values with `nmatlist2heatmaps()`
* Consider wrapper `nmathm2df()`

   * Convert coverage heatmap data to some type of `data.frame` suitable
   to save in Excel or tab-delimited format. Goal is to capture things
   like row clusters, row order, and row annotations, alongside
   the peak/row names.
   * Consider other helper functions to extract data or prepare useful
   data consistent with the heatmaps. For example, profile metaplots,
   so they could be plotted again separately.


# TODO 13nov2024

* Add GeoMx data processing functions.
* `design2colors()` - needs improvements - not nearly as automated as hoped

   * Fix error when supplying data that contains `NA` or empty values.
   * Fix warning:
   `"matrix(ncol=2, byrow=TRUE, add_hue), data length is not a multiple"`
   * Consider adding `default_colors` with common categories:
   Male/Female, TRUE/FALSE, Yes/No, positive/negative.
   Maybe `"Wildtype"`,`"WT"`?
   * Allow `sample_color_list` to inform "everything" when supplied.
   E.g. absence of `group_colnames`, `lightness_colnames`, `class_colnames`
   it should still be able to apply every column function.
   * `color_sub` - should accept `list` input, named by column,
   * `lightness_colnames` - should still use `color_sub` when matched?
   It should apply the lightness to the group color, if not the other colors.
   * When supplied group colors, consider defining class colors using group
   median color hue.
   * Consider option for columns to make their own categorical colors.
   * Consider matching cardinality of other columns to class, group, lightness,
   other columns with colors assigned.
   E.g. If "Genotype" and "Knockout" have values that are always 1:1 the same,
   it's good to share categorical colors.
   This could affect many columns, which is also fine.

* `rmd_tab_iterator()`

   * Consider `quiet=FALSE` default, which prints errors but continues,
   and is otherwise not verbose.


# TODO 08aug2024

* `rmd_tab_iterator()`

   * Consider option to print the variable name with the first tab:
   
      * Instead of: 
      ```
      [Expression], [Correlation]
      [Cohort], [Age]
      ```
      * It would be:
      ```
      [Data: Expression], [Correlation]
      [Centering: Cohort], [Age]
      ```


* `design2colors()`

   * `class_colnames` appears to be broken, it does not seem to assign
   class colors as expected.
   
      * It should use class to define range of color hues.
      * Group colors are assigned within each range of hues.
      * This part appears to be correct.
      * It should then assign class hue by the middle hue in the range,
      which does not appear to be correct.
      * It appears that class columns are assigned categorical colors
      as any other categorical column.
      * Also it defines `class_group_color` and `class_group_lightness_color`,
      but there should be `class_color` that corresponds to the one-or-more
      values in `class_colnames`.
   
   * Consider adding a function to plot the matrix separately, table format,
   as done by default in `design2colors()`. The `showColors()` and
   `print_color_list()` are not great when there is this much data.
   * It needs more intuitive "default behavior", for example define
   `class_colnames` then `group_colnames` and have it use the last value
   as `lightness_colnames`?
   
      * Can it use some logic to order factors by number of levels?
      Take only columns with at least 2+ replicated values, order from
      most to least replicated values to define:
      `class_colnames`, `group_colnames`, then `lightness_colnames`.

# TODO 10jul2024

* Consider new function `import_omics_data()`

   * Generic data matrix import routine.
   * Optional `curation_txt`, consistent with other import functions.
   * Column headers are assumed to be samples, only one row of headers
   is permitted (for now).
   * First column is expected to contain rownames.
   * There can be multiple row annotation columns, to be stored in `rowData()`.


# TODO 28jun2024

* Consider function to display a `data.frame` using `sample_color_list`
output from `design2colors()`.

   * It should essentially do the same thing as `design2colors()` except
   re-use the `sample_color_list`

# TODO 12jun2024

* Consider generic importer for tab-delimited, csv-delimited data
that produces a `SummarizedExperiment` object.

   * Bonus points for some mechanism to recognize multiple header rows,
   so they can be added to `colData()`.
   * Bonus points for option to specify whether rows or columns are
   measurements or samples, i.e. `rowData()` or `colData()`.

# TODO 03jun2024

* `nmatlist2heatmaps()`

   * Consider new argument `max_partition_rows` to limit the number
   of rows for any given partition. This option needs to be done here
   particularly so it can be used together with `k_clusters` which
   produces clusters whose sizes are not known beforehand.
   * Consider option to rotate color legend labels 90 degrees (vertical)
   so they take less width.

# TODO 30may2024

* Coverage heatmap automation functions. So far they include:

   * `make_cov_config()` - create `config_df` for matrix files
   * `auto_cov_heatmaps()` - given `config_df`, "draw the rest of the owl"
   (Do all the things required to make a consistent set of heatmaps.)
   * `rbind_cov_config()` - to combine multiple independent `config_df`
   * `process_cov_files()`
   
      * Confirm matrix files exist
      * When provided coverage (and not matrix), call `make_heatmap`
      to create matrix files
      * No doubt this is a stop-gap specific to certain systems, e.g.

         * call `system()` or `system2()`
         * call `sys::exec_wait()` or `sys::exec_background()`
         * call `processx::run()` - preferred for multiple concurrent processes


# TODO 20may2024

* `parse_ucsc_gokey()`

   * DONE. Consider mechanism to define custom default track settings.
   For example `get_track_defaults()` could be used to define
   settings in the package environment, which could then be edited.
   Similar to `igraph::add_shapes()` and `igraph:::.igraph.shapes`,
   an `environment` which is created during package loading.
   This `environment` can be adjusted.

* Add `nmat2coverage_matrix()`

   * Intended to save a matrix file in tab-delimited format.
   * Workflow would be to manipulate the matrix, possibly calculating
   sum, mean, subtraction, etc. Then save the result to a file for
   easier, more consistent re-use.

* `nmatlist2heatmaps()`

   * DONE. Debug why the `k_colors`, `color_sub` usage is not working properly
   with `row_split`.
   * Handle situation when input data is all zero. The detected range
   should be set to c(0, 1).
   * Handle empty matrix (entirely NULL) as if it were zero.

* `coverage_matrix2nmat()`

   * Consider handling missing file by creating a matrix with zeros?
   Use case is to intentionally create empty matrix to be filled with
   derived data, for example matrix 3 minus matrix 2; or mean of matrix 1
   through matrix 8.

# TODO 17may2024

* DONE. `nmatlist2heatmaps()`

   * DONE. Consider applying word-wrap to color legend titles, so the title
   does not become unreasonably long.

# TODO 16apr2024

* DONE. `rmd_tab_iterator()`

   * DONE. The `test=TRUE` functionality is not working properly for nested tabs.

# TODO 10apr2024

* `design2colors()`

   * Consider method to assign categorical colors to a column name,
   for example `Sample_ID="rainbowJam"` would create categorical colors
   for all values in this column, preventing it from assuming the
   same color as the sample group.
   * Consider allowing `color_sub` to be a `list` which may contain
   a `function` to use for color assignment of column values.

# TODO 01mar2024

* `get_salmon_meta()`

   * when a JSON file is missing, give an error with that information,
   for example when `meta_info.json` is absent, the error should
   say that, to help the user find and fix the problem.

# TODO 26jan2024

* `rmd_tab_iterator()`

   * Consider option to define the active tab for each layer of tabs,
   see https://bookdown.org/yihui/rmarkdown-cookbook/html-tabs.html
   `### Tab name {.active}`

# TODO 15dec2023

* `import_salmon_quant()`

   * Done. Debug and fix occasional error `"duplicate 'row.names' are not allowed"`

# TODO 06dec2023

* migrate `slicejam::import_featurecounts()` here

   * add argument `curation_txt` to populate `colData()` consistent with
   other uses in this package.

## `rmd_tab_iterator()`

* allow tabs to be hidden

   * DONE. design idea: optional argument to `base_fn` `test=TRUE` will return
   `logical` indicating whether to display the tab.
   * DONE. When `test` argument is defined for `base_fn()`, for example
   `base_fn <- function(..., test=TRUE){ 1 }`:
   
      1. First call: `base_fn(test=TRUE, ...)`
      2. If `FALSE` the tab is hidden and proceeds to the next iteration.
      3. If `TRUE` the tab is shown, then call: `base_fn(x, test=FALSE, ...)`
      4. If the return value is not `logical`, we assume the tab contents
      have already been displayed.
   
   * When `test` argument is not defined for `base_fn()`, for example
   `base_fn <- function(...){ 1 }`
   
      1. Tab is displayed, and a call is made to `base_fn(...)`

* catch errors in `base_fn` so the tabs will continue to iterate

   * consider returning `logical` for success (`TRUE`) or any error (`FALSE`),
   or `integer` number of errors, so `0` indicates zero errors, `4` indicates
   four total tabs caused an error.
   * consider `tryCatch()` around the `base_fn()` calls.

# TODO 09nov2023

* `rmd_tab_iterator()`

   * it has been useful to wrap `base_fn` inside `tryCatch()` which
   continues and prints the error, without crashing the RMarkdown.
   Consider how this may be added into this function as a default action.
   * Need option to hide a tab, most commonly when there is no suitable
   plot to be produced given the combination of parameters.
   
      * Design idea: return value `FALSE` indicates the tab should not
      be displayed? Unclear how the order of steps might permit
      hiding a tab...

* `curate_se_colData()`

   * Add an option to reorder `se` samples based upon the order found
   in the curation `df` data.

* `import_nanostring_csv()`

   * DONE. Update `import_nanostring_csv()` to subset by matched rows
   in `curation_txt`, and reorder samples to match `curation_txt`.

* `import_nanostring_rcc()`

   * add argument `curation_txt` to behave as `import_nanostring_csv()`

# TODO 27oct2023

* `nmatlist2heatmaps()`

   * **Migrate into coverjam**
   
      * replace these functions with dependency on Github `"jmw86069/coverjam"`

   * more examples of customizing font sizes
   * consider option for global adjustment of all font sizes together
   * consider `legend_base_nrow` logic that inspects the label width,
   to avoid multiple columns when labels are already very wide, it squishes
   the coverage heatmap size making the heatmaps too narrow.
   * investigate whether figure size can be calculated/predicted upfront,
   instead of having to calculate by `length(nmatlist)` and estimating the
   padding required for row annotations, heatmap gap spacing, and however
   large the color legends might be.

* new functions

   * `subset_nmatlist()` - convenience function to subset all matrices
   in the list. Intended for row subsetting `rownames(nmatlist[[1]])`
   but could also use nmat zoom options.

* consider `normalizedMatrixList` object class?

   * it would extend `EnrichedHeatmap::normalizedMatrix` and enable subsetting
   functions such as `nmatlist[x, ]` which thereby subsets rows in all matrices
   * it could become formal input to `nmatlist2heatmaps()`
   * it might be able to store associated heatmap parameters:

      * `signal_ceiling`
      * `ylim` used by `anno_enriched()`
      * `nmat_color` color function or color ramp
      * `label` used as title above the heatmap
   
   * it could simplify pre-processing, such as creating "diff" matrices,
   or any combination of matrices in the original `nmatlist`.

# TODO 23oct2023

* `nmatlist2heatmaps()`

   * DONE. consider storing important function arguments in the returned object
   
      * `k_clusters`, `k_method`
      * `rows`
      * `anno_df` (?) it could be a large data object, or `colnames(anno_df)`
      * `color_list` - the list of colors used for `anno_df`
      * `nmat_colors` (as the final color function for each heatmap)
      * `panel_groups` and the parameters for each panel group:
      
         * `ylims`
         * `signal_ceiling`

   * DONE. Consider figure caption with clustering information, similar to
   the caption used by  `multienrichjam::mem_gene_path_heatmap()`:
   
      * caption is stored as `attr(hm, "caption")`
      * the function that renders the caption is also stored:
      `attr(hm, "draw_caption")`
      * the caption is rendered using `grid::grid.textbox()`
      * suggested caption contents:
      
         * `"N rows displayed"` OR
         * `"N rows partitioned into M groups"`
         * `"k-means method='correlation' using heatmap X"` OR
         `"k-means method='correlation' using heatmaps X,Y,Z"`
         

# TODO 13oct2023

* `rmd_tab_iterator()`

   * Need way to "hide" tabs:
   
      * when no plot is produced
      
         * Can code determine whether a plot is produced? probably not easily.
      
      * certain combinations of tab values that should be hidden
      
         * optional accessory function that returns TRUE or FALSE, on whether
         to call `base_fn`.
      
      * driving example: MA-plots with `useRank=TRUE` are identical
      for raw and normalized data, should be shown once, for raw data,
      optionally for batch-adjusted data. Either need the MA-plot code to
      "hide" unnecessary plots, or provide criteria for which tabs to hide.

   * consider option to display heading label before each set of tabs

      * Previous markdown example does not include the tab type, just
      the tab value:
         ```
         # figures
         
         ## tab_list_1 value_1 {.tabset}
   
         ### tab_list_2 value_1
         
         {figure 1}
         
         ### tab_list_2 value_2
         
         {figure 2}
   
         ## tab_list_1 value_2 {.tabset}
   
         ### tab_list_2 value_1
         
         {figure 3}
         ```

      * New markdown example includes the tab type first, then the next
      layer uses each tab value:
      ```
      # figures
      
      ## tab_list_1 label {.tabset}
      
      ### tab_list_1 value_1
      
      #### tab_list_2 label {.tabset}
      
      ##### tab_list_2 value_1
      
      {figure 1}
      
      ##### tab_list_2 value_2
      
      {figure 2}

      ### tab_list_1 value_2

      #### tab_list_2 label {.tabset}
      
      ##### tab_list_2 value_1
      
      {figure 3}
      ```

# TODO 26sep2023

* migrate `slicejam::import_featurecounts()`

   * currently only imports the file format into `data.frame`
   * goal is to create fully described `SummarizedExperiment`
   using `rowRanges()` for genomic coordinates.
   * add argument `curation_txt`

* migrate `slicejam::fc_to_curation_txt()`

   * convert input data into a reasonably good `"curation.txt"`


# TODO 21sep2023

* DONE. Update `design2colors()` consistent with colorjam 0.0.26.900 changes.

# TODO 31jul2023

* `import_salmon_quant()`

   * be more tolerant of being given paths to `quant.sf` files which
   may be renamed to different filenames.
   * be more tolerant when GTF and tx2gene are not supplied, the import
   can still import transcript-level data without adding gene-level
   annotation and summaries.

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
