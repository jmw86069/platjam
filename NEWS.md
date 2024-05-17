# platjam 0.0.82.900

* Added `stringi` to package dependencies, for word wrap without incurring
heavy dependencies from `stringr` (which calls `stringi` anyway).

* `nmatlist2heatmaps()`

   * Added option to hide the caption, but `do_caption=TRUE` is default.
   * New arg `caption_fontsize=10` to control caption text size.
   * The caption includes `byCols` columns in the `anno_df` used to
   order the rows.
   * When data is transformed, the transform is no longer printed
   inside parentheses.
   * The color legend now uses `adjust_annotation_extension=TRUE` to
   ensure the legend does not overlap the annotation labels.
   * Color legend titles are word-wrapped at 15 characters, to handle
   very wide multi-word labels.
   * `raster_by_magick` by default checks that `magick` is available.

## new functions

* `apply_word_wrap()` - calls `stringi::stri_wrap()` and combines
vectors by `sep`, default `"\n"` for multi-line word wrap by default.
It improves `base::strwrap()` by also offering `whitespace_only=TRUE`
so that strings can be split at `"-"` and other light punctuation,
instead of requiring a space `" "`.


# platjam 0.0.81.900

## bug fixes

* `rmd_tab_iterator()`

   * Fixed bug when `base_fn` contained argument `test`, the environment
   was not properly defined before testing whether the tab should be
   visible, even though the environment was properly defined before
   calling the subsequent function `base_fn(..., test=FALSE)`.

# platjam 0.0.80.900

## bug fixes

* `import_nanostring_rcc()`

   * added prefix to call `jamba::pasteByRow()`
   * changed detection of positive/negative control probes to use `"^POS_"`
   and `"^NEG_"` instead of previous `"^POS"`, `"^NEG"`. Previous approach
   would mark `"POSTN"` a positive control, when it is an endogenous gene.

## changes to existing functions

* `import_nanostring_rcc()`

   * New argument `curation_txt` to include sample curation during import.
   
* `import_nanostring_rlf()`

   * changed default `plot_type="none"` so the default is not to create
   a color assignment plot.

# platjam 0.0.79.900

* `rmd_tab_iterator()` - major update

   * `base_fn()` can now contain optional argument `test=TRUE` to determine
   whether a particular tab should be displayed.
   
      * When `base_fn()` is called with `test=TRUE` it is expected to
      return `logical` indicating whether the tab should be displayed,
      and should provide no other output. It will be called again
      with `test=FALSE` if the tab should be displayed, in which case
      the function should generate output.
      * When `base_fn()` does not contain named argument `test`, the
      tab is always displayed, the function is called only once, and
      it should generate output as normal.

   * `base_fn` and `fn_list` functions are all wrapped inside `tryCatch()`
   which prints the error message and does not exit the iterator.
   * `fn_list` will likely be deprecated and retired in the near future,
   relying upon `base_fn` for all tasks.

* `import_salmon_quant()`

   * Fixed rare error `"duplicate 'row.names' are not allowed"` caused
   by `"gene_body"` transcripts (unspliced transcripts) which also
   contained the strand `"(+)"` or `"(-)"` and were being removed
   by the default argument: `trim_tx_from=c("[(][-+][)]")`.
   This situation occurs when one gene locus (defined by `"gene_name"`)
   is present on two strands (which probably should never happen,
   but in some Gencode files it does; perhaps by using `gene_name` and
   not `gene_id`). The argument default is now `trim_tx_from=NULL`
   which also makes sense that the default condition is not to change
   the rownames, nor the values in `txColname="transcript_id"`.
   * Fixed argument help docs, which incorrectly used `curate_tx_from`,
   `curate_tx_to` instead of `trim_tx_from`, `trim_tx_to`.

* `print_color_list()`

   * fixed error with function elements using the updated hex style
   from `circlize::colorRamp2()`.
   * added examples.


# platjam 0.0.78.900

* `import_nanostring_csv()` two changes when supplying `curation_txt`:
   
   1. The output `se` is subset to include only sample which were
   annotated in the `curation_txt` file.
   2. The order of samples now matches that with the `curation_txt` file.

# platjam 0.0.77.900

* `nmatlist2heatmaps()`

   * new argument `anno_row_gp` with default fontsize 14 for row
   marks when `anno_row_marks` is provided.
   * `anno_row_labels` now accepts one or more columns in `anno_df`
   to define labels, when `anno_row_marks` is provided. Not super convenient,
   because row labels are most likely not shown as `anno_df` columns.
   * `anno_row_labels` now tolerates when `names(anno_row_labels)` is not
   defined, and uses `anno_row_marks` by default.
   * apply arg `axis_name_gp` to `anno_enriched()` axis labels, for consistent
   font size for x-axis labels, and profile y-axis labels.
   * `axis_name_gp` default fontsize 10, instead of 8 previously.
   * `column_title_gp` default font size is 14, instead of 12 previously.
   * `legend_base_nrow=12` is default, mainly to provide two columns for
   human chromosomes, but one column for most other scenarios.
   * new argument `title_gp` with default fontsize 14 for the overall title.


# platjam 0.0.76.950

## changes to existing functions

* `zoom_nmat()` now accepts `""` or `"none"` to avoid transformation
* `nmatlist2heatmaps()`

   * new default `main_heatmap=NULL` causes all heatmaps to be used for
   row ordering. The enriched score is calculated for each heatmap,
   summed by row, then used to produce decreasing order.
   * new argument `transform_label` to customize transformation labels
   included below each coverage heatmap title. By default it uses
   the `character` string from `transform` for anything except `"none"`.
   * new argument `padding` to add whitespace padding around the heatmap set.
   * `axis_name` now accepts single `character` input to customize the
   center label.
   * A default caption is added to the bottom-right corner of the output,
   which describes basic information about the figure: number of heatmaps,
   number of rows, whether k-means was used and which method, the
   total number of row partitions, and which heatmaps were used for
   row ordering and k-means clustering.
   * returned data include `fn_params` with a `list` of important parameters
   associated with the heatmaps. Notably some values are calculated
   on the fly, like when `panel_groups` is defined: `ylims`, `signal_ceiling`,
   `nmat_colors`, etc.
   * default color ramp is `"Reds"`
   * `main_heatmap` can now accept multiple values, causing `row_order`
   by default to order using the sum `EnrichedHeatmap::enriched_score()`
   across each heatmap, before being ordered as usual.
   * `row_order` will by default substitute `NA` values with `0`, in rare
   cases that `NA` values may exist, maybe caused during transformation.
   * argument `k_method="correlation"` is a new default, which falls back
   to `"euclidean"` if the `amap` package is not installed.
   * new argument `min_rows_per_k=100` requires at least 100 rows per k-means
   cluster, to protect from clustering with a very small number of rows.
   The argument was motivated by using `partition` and `k_clusters` together,
   for which the number of rows in each partition is not (easily) known
   upfront. However, it is also useful to set `k_clusters` and be assured
   that at least `min_rows_per_k` on average are available during clustering.
   * Added new example demonstrating use of `partition` and `k_clusters`.
   * Added and updated function documentation for clarity.
   * silenced verbose output `"Preparing ComplexHeatmap::draw()"`
   * added more error-checking and specific error messages when matching
   across the optional inputs, all of which must match when provided:
   
      * `rows`
      * `rownames(nmat)`
      * `rownames(anno_df)`
      * `names(partition)`

   * `partition` is handled internally as a `factor` so the level order
   is maintained. If supplied as `factor` the levels are honored as given,
   even when combined with `k_clusters` as below.
   * `kmeans()` warnings are suppressed, otherwise they interfere with
   RMarkdown output.
   * `k_clusters` and `partition` work together
   
      * each partition is individually k-means clustered
      * `k_clusters` can be a vector, applied to each `partition` in order,
      or directed by `names(k_clusters)` if all partitions have an associated
      value in `names(k_clusters)`. Otherwise `k_clusters` is recycled to
      the number of partitions, then applied in order.
      * For each partition, `k` requires at least 10 rows per `k` rounded up.
      A partition with 10 rows can only have `k=1`, and partition with 11 rows
      can have `k=2`, etc. This step protects from k-means clustering small
      partitions.
      * partition colors are assigned to `partition` first, then colors are
      split by `jamba::color2gradient()` across k-means clusters within each
      partition.
   
   * `k_method` now also offers `"spearman"` for rank-based correlation metrics,
   although this method has not been widely tested on genome coverage data.

# platjam 0.0.75.900

## new functions

* `rmd_tab_iterator()`

   * Intended to help RMarkdown documents generate tabbed output,
   even allowing multiple layers of tabs,
   usually to configure the same visualization with a few different options.
   * Provide a list of vectors to use for tabs.
   * It iterates each layer of tabs, defines variables in the appropriate
   environment to be "seen" by subsequent layers of tabs.
   * The base layer calls a function, usually producing a figure or table.
   * The process clean up R code, preventing multiple layers of `for` loops
   or `lapply()` calls.
   * This function may be moved into `jamba` eventually, but for now is
   being tested and may be updated rapidly based upon real use cases.


# platjam 0.0.74.900

* `design2colors()`

   * simplified how class and group colors are assigned overall
   * simplified the call to `colorjam::rainbowJam()`, no longer
   assigning color hue (then applying color warping in colorjam),
   instead using `colorjam::rainbowJam()` since it has additional
   logic for low `n` color hue assignments. Pushes all aesthetic
   choices to `rainbowJam()`.
   * Changed how `phase` is applied, so the ordered phase is applied only
   to non-padded colors. When `class_pad` or `end_hue_pad` insert blank
   colors to absorb the corresponding color hue, the phase is not consumed,
   since phase is intended to be applied across the final set of colors.
   Previously, the phase would be out of phase with the final colors,
   so to speak, making the output colors appear weird.
   * Fixed application of `end_hue_pad` and `class_pad` together, previously
   it only applied `class_pad`.
   * `class_pad` is applied partially at the start and end of the color wheel,
   splitting the padding before the first and after the last color hue.
   * `color_sub` can substitute for `class` or `group` colors without
   requiring all values in the column be matched. Useful for altering one
   color in a set.
   * Improved color saturation for `class` colors created when blending
   corresponding group colors. The process uses
   `colorjam::vibrant_color_by_hue()` to determine the most vibrant color
   for the same hue, then blends this color with the original class color.

# platjam 0.0.73.900

* dependency version bump for colorjam (>= 0.0.26.900)

## changes to existing functions

* `quick_complement_color()` was updated to match recent updates to colorjam
version 0.0.26.900. This change also affected `design2colors()`
* `design2colors()`

   * `preset="dichromat2"` new default, so colors begin with yellow,
   more conducive for applying to control groups than starting with red.

# platjam 0.0.72.900

## migration of coverage heatmap functions to coverjam

The functions related to coverage heatmaps are being migrated into
`"jmw86069/coverjam"` Github-hosted R package, for cleaner maintenance
of those functions.

# platjam 0.0.71.900

## changes to existing functions

* `import_metabolomics_niehs()`

   * added optional step to import `"compounds_pos.txt"` and
   `"compounds_neg.txt"` files, which contains important annotations
   for which measurements are imputed by the upstream software, and
   using which impute method.
   * When the files above are not present, and the file
   `"1_DataProcessed.zip"` exists, it is used to extract the raw files
   as named above.
   * When no matching files are detected, this import step is skipped.

## new functions

* `convert_imputed_assays_to_na()`

   * utility function used to apply filtering based upon the impute flags,
   thereby converting imputed data matching the relevant flags to `NA`,
   and storing in a new `assays()` entry for downstream analysis.

* `process_metab_compounds_file()`

   * Internal function used to import full compounds data supplied
   as a `data.frame` after loading from the `"1_DataProcessed.zip"` file.
   * The function will merge data into a `se_list` object if supplied,
   otherwise it will return the full `se` object including each assay
   name imported from the file.


# platjam 0.0.70.900

## bug fixes

* Added missing package prefix for several calls to `jamba::nameVector()`,
or `jamba::printDebug()`. Affected functions:

   * `import_nanostring_rcc()`
   * `coverage_matrix2nmat()`
   * `frequency_matrix2nmat()`
   * `nmatlist2heatmaps()`

# platjam 0.0.69.900

* `design2colors()`

   * input recognizes `DataFrame` which is commonly produced from
   `SummarizedExperiment::colData()`.
   * input recognizes `tbl_df` (tibble) and `matrix` by coercing to
   `data.frame`.
   * column sort can be reversed with "-" prefix. For example
   `"-Treatment"` will reverse the sort order in column `"Treatment"`.
   * allow empty `group_colnames` so colors are assigned in all columns
   * tolerate various combinations of `group_colnames`, `class_colnames`
   with non-ideal cardinality, by using different workarounds strategies
   * new arguments `Crange`,`Lrange` which are passed to
   `colorjam::rainbowJam()` to set a less saturated color palette by default.

# platjam 0.0.68.900

## bug fixes

* `import_metabolomics_niehs()`

   * Fixed bug when importing files that may contain duplicated
   rows in metadata, and in each data matrix file.
   The issue appears to be pure duplication of rows, so the workaround
   is to print a warning, and retain non-duplicated records
   per SampleID value.

# platjam 0.0.67.900

## new functions

* `import_metabolomics_niehs()`

   * imports metabolomics data specifically in file formats used
   by the NIEHS Metabolomics Core Facility. The output is a `list`
   of `SummarizedExperiment` objects, typically including positive
   and negative ionization data.

* `curate_se_colData()`

   * wrapper function to `curate_to_df_by_pattern()` specifically
   for `SummarizedExperiment::colData()`, updating data in place.
   * optionally this step can subset columns for those that
   match curated patterns.

# platjam 0.0.66.900

## bug fixes

* `import_salmon_quant()` was using the older version of
`curate_to_df_by_pattern()` supplied by `slicejam` instead of
the slightly updated version in `platjam` (this package).
Not such a big issue except that `slicejam` is not yet public.

# platjam 0.0.65.900

## changes

* added genejam Github information to DESCRIPTION, previously this
package was not properly referenced as part of Github and would
cause the package dependency to fail.

# platjam 0.0.64.900

## changes to existing functions

* `import_nanostring_csv()`

   * new argument `curation_txt` for optional sample annotations during
   import; behaves exactly as with `import_proteomics_mascot()`.
   * TODO: will add this argument to `import_nanostring_rcc()`
   * TODO: will create generic function for curation into SE `colData()`
   * fixed bug not prefixing `rowData()`

* `import_proteomics_mascot()`

   * updated to fix error when colnames are changed by `curation_txt`


# platjam 0.0.63.900

* added testthis unit tests for some functions

## changes to existing functions

* `assign_numeric_colors()`

   * new argument `color_max` to define a fixed upper `numeric` threshold
   for the maximum color in the color gradient.

* `design2colors()`

   * Changes to enable `color_sub` to accept color ramp names in addition
   to accepting single colors.
   * Slightly cleaned up some logic during color assignment.
   * New argument `color_sub_max` intended to be passed to
   `assign_numeric_colors()` argument `color_max` for particular
   column names.

# platjam 0.0.62.900

## changes to existing functions

* `import_proteomics_mascot()`

   * Added more robust detection of changes in column order between
   the input Excel file and the curation_txt data.
   * The output column order now matches the order of the curation_txt
   file, which should be helpful when defining the control group
   in statistical contrasts.

* `curate_to_df_by_pattern()`

   * new argument `order_priority="df"` to indicate the order of rows
   in the `data.frame` returned will match the curation `df`.
   * option `order_priority="x"` will order rows in the `data.frame`
   returned by the order of entries in `x`.

* `import_lipotype_csv()`

   * minor changes to ensure the order of curated rows and sample identifiers
   were matched when the curation data did not contain `"Filename"` column,
   by default uses the last column returned by `curate_to_df_by_pattern()`.

# platjam 0.0.61.900

## new functions

Proteomics data is sometimes provided in separate results files,
which can later be merged into one larger format for downstream
analysis.

* `merge_proteomics_se()`

   * This function provides a specific method of merging two
   `SummarizedExperiment` objects, with heuristics appropriate
   for proteomics expression data.

# platjam 0.0.60.900

* `jamba` dependency was bumped to 0.0.89.900
* added `SummarizedExperiment` and `XML` to dependencies

## updates to existing functions

* `import_nanostring_csv()` and `import_nanostring_rcc()`

   * both functions now return data transformed with `log2(1 + x)`
   * CSV import will name the assay `"norm"` if the CSV file contains
   that substring, otherwise defaults to `"raw"`, consistent with
   RCC import.
   * The `assay_name` can be defined upfront for both functions.

## new functions in development

* `validate_heatmap_params()`

   * This function begins a process of making the coverage heatmaps
   into a command-line tool. Work in progress.

# platjam 0.0.59.900

## updates to existing functions

* `import_salmon_quant()`

   * new arguments `curate_tx_from`, `curate_tx_to` to allow custom
   curation of transcript identifiers, specifically for transcripts
   where the FASTA workflow introduces `"(-)"` and `"(+)"` to the
   identifiers. This situation occurred when adding gene body transcripts
   to the Thuman T2Tv2.0 genome GTF entries.

* `design2colors()`

   * argument `color_sub` will convert color names to hex format, for
   compatibility with some downstream uses, for example `knitr::kable()`
   and `jamba::kable_coloring()` when used in an Rmarkdown document.


# platjam 0.0.58.900

## bug fixes

* `import_salmon_quant()` was errantly checking the class of `tx2gene`,
causing it to ignore user-supplied `tx2gene`. This issue affected import
when `gtf` was not supplied, which was much more common.

## changes to existing functions

* `import_salmon_quant()`

   * the list of assays now use names directly from `tximport()`, specifically
   `c("counts", "abundance", "length")`, and not `c("counts", "tpm")`.
   * `TxSE` now contains assayName `"length"`, and metadata item
   `"countsFromAbundance"`, so the output can be used with
   `tximport::summarizeToGene()` as a second step outside this import
   process.
   * new argument option `import_types="gene_body"` designed specifically
   when the GTF or `tx2gene` contain entries representing unspliced
   multi-exon genes, in addition to spliced proper transcripts. When present:
   
      * `import_types="gene"` will summarize transcripts to gene level,
      excluding entries annotated `"gene_body"`, returned as `GeneSE`.
      * `import_types="gene_body"` will summarize transcripts to gene level,
      including entries annotated `"gene_body"`, returned as `GeneBodySE`.
      * `import_types="gene_tx"` will summarize transcripts to gene level,
      and separately represents `"gene_body"` entries,
      returned as `GeneTxSE`.

# platjam 0.0.57.900

## new functions

* `import_nanostring_csv()` is a basic import function for the common
CSV style export from Nanostring software. This version has minimal
extra capabilities, and produces a simple `SummarizedExperiment` object.


# platjam 0.0.56.900

## changes to existing functions

* `import_proteomics_PD()` was updated to fix an errant check of
the protein colnames. This check occurs to verify that the accession
colnames did not contain delimited values, which causes the colnames
to be split into individual values by `jamba::makeNames()` with suffix.
When that happens, the original colname is added before the split columns
so the original data can be verified.

## new functions

* `import_proteomics_mascot()` is a version 1 importer for MASCOT and
similar data, which expects protein-level data with `totalIntensity`
and `numSpectra` columns for each sample.


# platjam 0.0.55.900

## bug fixes

* `import_proteomics_PD()` was updated to remove dependency on `slicejam`
R package which is currently private during internal development.
The only dependent function `curate_to_df_by_pattern()` was moved into
this package.

# platjam 0.0.54.900

## bug fixes

* `design2colors()` caused an error when there were no leftover
categorical colors to be assigned after the initial group/class/lightness,
then user-defined colors by colname. This bug has been fixed.

# platjam 0.0.53.900

## new functions

* `import_nanostring_rlf()` is an import function to load Nanostring
codeset data provided in the form of an RLF file.

## changes to existing functions

* `design2colors()` now calls `colorjam::blend_colors()` to define
colors for each class. The method is still in development, so the
specific methods may change. The output is intended to represent
the group colors contained in each class, we'll see how well it fares.

## bug fixes

* `design2colors()` was updated to convert columns with class `"table"`
into `"integer"`, otherwise weird things happened during coersion,
causing an error.


# platjam 0.0.52.900

## new functions

* `mean_hue()` is a simple function to return a mean hue in degrees.
* `print_color_list()` is a convenience function to convert a color list
that may contain color vectors, or color functions, into a list of
named color vectors. It will print to the console using those colors.

## changes to existing functions

`design2colors()` was updated to to handle an edge condition.

* Refactoring the workflow:

   * define colors by design:
   
      * define colors to `group_colnames`
      optionally partitioned by `class_colname`
      * define gradient colors to `lightness_colnames`
      * define class colors using the mean hue for groups contained in each class
   
   * Apply design colors to all matching columns:
   
      * Test all columns for compatible cardinality with group, lightness.
      * Note: column names in `color_sub` will be ignored at this step,
      instead those columns use `color_sub` to generate a gradient.
      * Assign colors to each compatible column.
      * Note this step will include extra annotation columns
      that have consistent cardinality.
      * (This might be the coolest step tbh, since it figures out which columns
      match the design.)
   
   * Expand `color_sub` with newly defined colors.
   
      * Intent is to assign consistent categorical colors to values
      present in other columns.
   
   * For columns without assigned colors:
   
      * Apply `color_sub` to the column when it matches the column name.
      * Numeric columns produce a color function.
      * Character columns produce a named color vector.

   * Assemble all remaining unassigned values:
      
      * Numeric columns use the column name as its "value", to define
      one color to the column header.
      * Remaining `character` values defined in `color_sub` are assigned.
      * Remaining unassigned `character` values, and `character` column names
      of numeric columns are assigned categorical colors.
   
   * Apply new categorical colors

      * For `numeric` columns, generate color function using the categorical
      color assigned to that column name.
      * For `character` columns, generate a color vector.

* Note that the mechanism used may change, I am still testing how
the assumptions work with real world data.

* New behavior: When `color_sub` is assigned to a value in `colnames(x)`,
the cardinality is not checked for consistency with group color assignments,
and instead will use the color assigned for that column.
* New argument `force_consistent_colors`:

   * During color assignment, sometimes a column value is assigned one
   color in a column, and a different color in another column. It usually
   happens when `color_sub` defines a color for a column, where that
   column unintentionally shares a value with another column that
   gets colors assigned earlier in the chain.
   * New argument `force_consistent_colors=TRUE` forces values to be
   assigned only one color. It seems obvious to use this option, but why
   not always?
   * For purely text strings, like "Treated" and "Untreated" one would
   think any time "Treated" appears in a graphic, it should use the same
   color. Yes. However, sometimes the number "2" appears in a column,
   as character or integer, it doesn't matter for this purpose. I have
   seen something like "sample_num" with values "1", "2", "3", etc.
   Their cardinality matches 1:1 with class/group/lightness color assignments
   defined by this function, so these values are usually also assigned
   to the group colors. However, these values may appear in another column
   defined by `color_sub`, intended to apply a gradient of that color to
   values in the column. What if that column contained values
   `c(2, 20, 500, 1000)`? Only these values would be colorized: `c(20, 50, 1000)`,
   while `2` would be assigned the color from group colors. It looks weird.
   * There are workarounds:
   
      * pre-filter columns which might have this effect. In the case above,
      remove `"sample_num"` before calling `design2colors()` so that column
      is not assigned a color.
      * associate all numeric columns with `color_sub`
   
   * The end result, which is intended, is to allow columns to have
   independent color assignments when values in those columns are not
   related. In the case above, the number `"2"` means something different
   in one column than another column.

# platjam 0.0.51.900

## changes to existing functions

* `design2colors()` was updated to adjust the low end of color gradient.

   * specifically when the low-high numeric range is above zero, for example
   `c(20, 50)`. The adjustment lowers the low end of color scale 5% the diff,
   so the low actual color will not be white.
   * Also the function was using `pretty()` to generate reasonable numeric
   steps across the numeric range, and was using those steps to define the
   color gradient. This method had weird effects, partly due to `pretty()`
   making choices about the best visible number for display, instead of
   the point range. Instead, the color gradient is applied to the calculated
   number range, then that function is used to generate a color function
   across the `pretty()` scale. I'm not sure it is ideal when `pretty()`
   does not color the numeric range, but it seems reasonable so far.


# platjam 0.0.50.900

## changes to existing functions

* `import_proteomics_pd()` was updated to handle absence of `curation_txt`
where sample names and labels may or may not be renamed for each
abundance matrix.

   * Specifically `label_colname` is detected as `"Label"`, any column
   containing `"Label"`, then `"^Input$"`. Note `"Input"` is the default
   first colname produced when there is no `curation_txt` file, and
   should represent a unique column identifier.

* `design2colors()` was updated:

   * `class_colnames` is required to have 1-to-many cardinality
   with `group_colnames`, or throws an error. As a result, class
   is guaranteed to be associated only with one set of columns.
   * class is assigned a color hue using the mean color hue of
   the groups it contains.
   * `color_sub` will apply a gradient color to `numeric` column values.
   * TODO: add some ability to pre-define class color with `color_sub`
   * TODO: add some ability to define numeric colors or color function
   for columns that contain `numeric` data. They are currently
   handled as categorical values, which is not ideal.

# platjam 0.0.49.900

## changes to existing functions

* `save_salmon_qc_xlsx()` new argument `verbose=FALSE`, thus changing
the default which included some verbose output by default.
* `df_to_numcolors()` new argument `trimRamp` passed to
`jamba::getColorRamp()` to trim edge colors from linear color gradient,
thus making the colors less dramatic.

# platjam 0.0.48.900

## new functions

* `design2layout()` is experimental, and is intended to determine
appropriate plot panel layout dimensions, given experimental groups
and batch information. It should help align panels within
batch, then within group, into columns, so plots can be visually
compared.
* `cardinality()` takes two vectors and determines the data cardinality,
in terms of "one-to-one", "one-to-many", "many-to-many", "many-to-one".
This function may help detect when experimental factors are compatible
with color assignment by group or other factors.
* `df_to_numcolors()` is an experimental function only used to
create a color substitution vector for use with `kable_coloring()`.

# platjam 0.0.47.900

## new function `import_lipotype()`

* Intended to import CSV formats provided by LipoType lipidomics data
* Returns `SummarizedExperiment` object, with rowData and colData
assigned as appropriate.
* Next iteration will parse optional `curation_txt` to assign
sample annotation for experimental design, factors, etc.


## new function `design2colors()` (in development)

See examples `? design2colors()` for usage.

* Intended as an extension to `colorjam::group2colors()`
that also utilizes subgroup (smaller subsets of groups), and 
group class (which includes multiple groups together).
* If this function works well through early testing, it likely moves
into `colorjam` for broader re-use.
* I find myself assigning categorical colors by group, then using
`jamba::color2gradient()` to split colors into light-to-dark gradient:

   * time course, 0-hr, 2-hr, 4-hr, 8-hr (light-to-dark)
   * treatment: untreated, treated (light-to-dark)
   * for each group: batch A, batch B (light-to-dark)

* I find myself assigning these categorical colors to each SampleID,
or to subgroups, or other factors, so the colors are consistent when
used in tables, and `ComplexHeatmap` output.

   * assign colors to all annotation columns with appropriate cardinality

* Finally, other factors need colors, they should differ from group colors

   * for each individual factor, assign unique categorical colors
   * re-use colors when a factor level already has an assigned color

* Default output:

   * colors assigned for every column, and every column value of the
   input `data.frame` argument `x`.
   * Plots a graphical table of sample annotations showing color assignments
   using `jamba::imageByColors()`.
   * `list` of colors, named by `colnames(x)`, suitable for use with
   `ComplexHeatmap::HeatmapAnnotation()`.
   * optionally a color vector, named by factor level.
   * optionally returns a `data.frame` with colors, in same order as input `x`,
   suitable for use directly in `jamba::imageByColors()`.
   
* Relevant arguments

   * `group_colnames`: defines categorical colors, one per group
   * `lightness_colnames`: assigns light-to-dark gradient to categorical colors
   * `class_colnames`: arranges groups into classes, so groups with same class
   are assigned similar color hues, with spacing between classes.
   * `preset`: defines the color wheel used by `colorjam::rainbowJam()`:
   
      * `dichromat`: dichromat-friendly color wheel, largely removes green
      * `ryb1`: red-yellow-blue (common painter's color wheel for color mixing)
      * `rgb`: red-green-blue (default R, shows how much of the wheel is green)


# platjam 0.0.46.900

## new function

* `import_salmon_quant()`

   * does the job of `tximport::tximport()`, then associates
   sample annotations via `curation_txt` to produce
   `SummarizedExperiment` objects at transcript and/or gene
   level.

# platjam 0.0.45.900

## updates to existing functions

* `get_salmon_meta()`

   * now parses the fragment length distribution if
   the file exists. It returns mean, median, mode, and 25/75 percentiles.
   * Also the `names(metafile)` are maintaines as rownames in the output.

## new functions

* `parse_salmon_flenfile()` parses the Salmon fragment length file
usually named `"flenDist.txt"`, and return summary values: mean,
median, mode, 25/75 percentiles.


# platjam 0.0.44.900

## updates to existing functions

* UCSC track hub function `parse_ucsc_gokey()` was updated because
apparently `glue::glue()` added a new argument `.trim` to trim
whitespace from every output line AND made `.trim=TRUE` the
default value. Because why be consistent with past versions of `glue`.
Thank you again to all the R package authors that have chosen not
to change the output of the primary function in their package
by default. :) Not shown: All the updates I did not have
to do because other core functions keep working as normal.

# platjam 0.0.43.900

## updates to existing functions

* `nmatlist2heatmaps()` was updated to handle transformation of
data during k-means clustering. Previously this function was not
correctly defined. 

# platjam 0.0.42.900

## updates to existing functions

* `import_proteomics_PD()` gained some options to handle PD data content,
such as referring Protein gene annotations into the peptide data during
import. Apparently sometimes they only annotate the protein level, so
this step recovers some annotations otherwise lost.
Also sometimes peptide data is split across multiple rows when the
peptide may appear in multiple positions in the same protein.


# platjam 0.0.41.900

## new functions

* `import_proteomics_PD()` will import proteomics abundance data
from Proteomics Discoverer software, focused on the Excel format
in use in 2022 that includes protein abundance data, with nested
sub-tables of peptide abundance data with different columns.
The function will parse each type of row and return two
`SummarizedExperiment` objects, one for protein data, one for
peptide data. This function calls an internal helper function
`convert_PD_df_to_SE()` which performs the work of parsing
gene annotations, peptides, post-translational modifications,
sample annotations, and optionally reads `"curation.txt"` to
define specific sample group annotations for downstream
analysis.


# platjam 0.0.40.900

## changes to existing functions

* `parse_ucsc_gokey()` was updated to allow additional parameters to be
defined through `...` dot list. These values will override all values
defined in the track lines.
   * Examples have been added to demonstrate minimal use of the function,
   showing composite track views, and multiWig overlay tracks.
   * new experimental argument `multiwig_concat_header=TRUE` allows
   adjusting the format of multiWig superTrack names, being tested
   inside the UCSC genome browser.


# platjam 0.0.39.900

## changes to `nmatlist2heatmaps()`

* The `ComplexHeatmap::draw()` function by default includes the
argument `adjust_annotation_extension=TRUE` described here
https://jokergoo.github.io/ComplexHeatmap-reference/book/a-list-of-heatmaps.html#adjust-blank-space-caused-by-annotations

   * The problem was that y-axis labels on the right side of profile plots
   overlapped the color legends on the right side.
   * Alternatively, profile axis labels may be positioned on the left side,
   which avoids having to adjust the legend on the right side.

* New argument `top_axis_side` with values "right", "left", "both",
"all", "none"
* New argument `k_heatmap` which defines which heatmap or heatmaps are
used for k-means clustering. Key new feature, k-means clustering can
now involve multiple heatmaps, which may be useful to find patterns
across multiple signals. This feature is being tested currently.

* `save_salmon_qc_xlsx()` was updated to include jamba:: function
prefixes as needed.

## functions removed

* `set_xlsx_colwidths()` and `set_xlsx_rowheights()` were removed,
since they have been moved into package `"jamba"`.
* `twostep_gradient()`, `showDichromat()`, and `make_jam_divergent()`
were moved into package `"colorjam"`.
* `jam_divergent` and `jam_linear` were moved into package `"colorjam"`.


# platjam 0.0.38.900

## changes to `nmatlist2heatmaps()`

A series of changes was made to `nmatlist2heatmaps()`:

* New argument `top_annotation` to customize the annotation appearing
above each heatmap. It takes several forms:

   1. `TRUE` uses the default `ComplexHeatmap::HeatmapAnnotation(EnrichedHeatmap::anno_enriched())`
   with some platjam-specific customizations, mainly to match the profile line
   colors to the `partition` and `k_colors` arguments for consistency.
   2. `FALSE` hides the annotation
   3. `HeatmapAnnotation` uses the function as provided
   4. `list` of any combination of the above, recycled to `length(nmatlist)`.
   The `list` input can be used for specific plots above each matrix,
   which can also specify the data used for those plots. Note option
   3 above uses the closed form `EnrichedHeatmap::anno_enriched()` in order
   to be applied to the appropriate data matrix for each heatmap.

* New argument `top_anno_height` which is a `grid::unit` object specifying
the height of the `top_annotation` - only used when `top_annotation` is
not supplied in custom form.
* Help text for argument `signal_ceiling` includes reference to the
child function `get_nmat_ceiling()` which controls the logic:

   * `signal_ceiling` above `1` applies the strict numeric value
   * `signal_ceiling` between `0` and `1` applies a `quantile(x, probs=signal_ceiling)`
   * `signal_ceiling=NULL` uses the maximum absolute value

* Argument `row_order` is intended to work alongside `byCols`, so the
sort of `anno_df` will also inlude sort for row order within each subset
of `anno_df` as applicable.
* New arguments `raster_quality` and `raster_by_magick` passed to
`ComplexHeatmap::Heatmap()`.
* New arguments `upstream_length`, `downstream_length` used for optional
zoom of the coordinate range being displayed. These arguments are
passed to `zoom_nmatlist()` if either is not `NULL`.

* New behavior when `panel_groups` are supplied:

   * `heatmap_annotation_param` will automatically hide duplicate color legends,
   and will use the `panel_groups` value for the color legend title.
   * `top_annotation` default axis will be hidden except for the last
   in each set of adjacent duplicated panel_groups. For example,
   `panel_groups=c("A", "A", "B", "B", "A", "A")` will display y-axis for
   entries `c(2, 4, 6)`.
   * In future, heatmap panels may have slightly wider gap between
   different panel_groups values.

## new functions

* `zoom_nmatlist()` and `zoom_nmat()` are intended to help zoom into
a particular range of coverage matrix data, where the stored coverage
is wider than desired for the heatmap figure. They filter the columns
in each `normalizedMatrix` and adjust the attributes stored in:
   * `attr(nmat, "upstream_index")` - the column index positions upstream the target region
   * `attr(nmat, "downstream_index")` - the column index positions downstream the target region
   * `attr(nmat, "target_index")` - the column index positions representing the target region
   * `attr(nmat, "extend")` - the genomic distance upstream and downstream the target region



# platjam 0.0.37.900

## changes to existing functions

* `get_track_defaults()` was updated to adjust the default bigBed
track settings, so the default display is `"pack"` but can be adjusted
as needed.

# platjam 0.0.36.900

## changes to existing functions

* `parse_ucsc_gokey()` handles bigBed tracks with specific composite
track templates and relevant default values.

# platjam 0.0.35.900

## changes to existing functions

* `parse_ucsc_gokey()` was updated to improve the method of
naming supergroup and group tracks, so the supergroup
name will persist for every group it should contain.
* `nmatlist2heatmaps()` new argument `trim_legend_title` will
trim heatmap legend title to remove any text following
a newline. For example, the heatmap label may contain multiple
lines of text which are useful to the heatmap, but less useful
to the color legend. In that case, only the first line of text
is used as the heatmap legend title.


## new functions

* `nmathm_row_order()` returns the rowname order of the heatmap
produced by `nmatlist2heatmaps()`. When rows are partitioned,
the output is a `list` of rowname vectors.


# platjam 0.0.34.900

## new additions

* `jam_linear` and `jam_divergent` are new color gradients
intended to be visibly distinct, and color-blind-friendly.
The divergent scales in particular are designed so that
the two directions are visibly distinct even under three
different color blindness simulations from
`dichromat::dichromat()`.

## changes to existing functions

* `nmatlist2heatmaps()` changes that forced heatmap legend direction
to be `"vertical" were reverted, since Dr. Gu updated
`ComplexHeatmap` to fix the issue with `grid` on R version 3.6.*.
* `nmatlist2heatmaps()`:

   * new argument `pos_line` to control display of lines around the target region
   * return values include `HM_drawn` which is the object returned
   after `ComplexHeatmap::draw()`, and which contains all the column and
   row orders for each heatmap.

# platjam 0.0.33.900

## bug fixes

* `nmatlist2heatmaps()` encountered an error when displaying
heatmap legend with `direction="horizontal"` which appears
to be specific to R version R-3.6.1, and is not present
in R-4.0.0 and above. Specifically I believe this error
is from `grid` package version `3.6.1` and how the function
`unit.arithmetic()` is called by ComplexHeatmap. Nonetheless,
the workaround is to force `direction="vertical"` when
`grid` package version is below `4.0.0`.


# platjam 0.0.32.900

## changes to existing functions

* `nmatlist2heatmaps()` new argument `show_heatmap_legend` is
recycled to the `length(nmatlist)` to allow hiding legends
for certain heatmaps. This update is intended when heatmaps
are grouped by `panel_groups` which allows multiple heatmaps
to have the same shared numeric color mapping.

# platjam 0.0.31.900

## updates

* `nmatlist2heatmaps()` updated so `ylims` will take
priority over `panel_groups=TRUE` which defines its own
ylims ranges.


# platjam 0.0.28.900

## bug fixes

* `nmatlist2heatmaps()` bug with no `k_clusters`,
no `partition`, and no `k_colors` would fail during
profile plot because it had NULL color.

# platjam 0.0.28.900

## bug fixes

* `nmatlist2heatmaps()` does more work to validate `rows`
with respect to rownames of `nmatlist`, `rownames(anno_df)`,
and `names(partition)` to ensure everything is properly
synchronized. The smallest common set of `rows` will be
used.

## changes to existing functions

* `nmatlist2heatmaps()` now allows `partition` and `k_clusters`
to be combined, so rows can be both partitioned, and kmeans
clustered. Labels are created `"partition - cluster_number"`
and are ordered by partition then kmeans cluster number. So
if `partition` is a factor with ordered levels, that level
will be preserved in the output.
* `nmatlist2heatmaps()` new argument `color_sub` which allows
defined categorical colors by sending a named vector. Any column
in `anno_df` all of whose values are in `names(color_sub)` will
be colored accordingly.
* `nmatlist2heatmaps()` will use `k_subset` to re-order
partitions, in addition to plotting only the subset of
row partitions. The code was change to handle the
sequence of operations: kmeans, partitioning, k_subset,
then anno_df, marked rows.
* `nmatlist2heatmaps()` changed the default `ht_gap` down to 3 mm.

# platjam 0.0.27.900

## changes to existing functions

* Numerous small bug fixes in `nmatlist2heatmaps()`, mostly in
handling annotation categorical color legends.

# platjam 0.0.26.900

## changes to existing functions

* `nmatlist2heatmaps()` new default `transform="none"`, previously
was `jamba::log2signed()` which enforced log2 transformation.
Frankly not sure which is best, might depend upon the range of
data, which depends upon the type of coverage signal being used.
* `nmatlist2heatmaps()` now treats annotation columns in `anno_df`
as bi-directional when the max value is <= 50, which helps in cases
that have values `c(-1, 0, 1)` and sometimes only `c(0, 1)`.
The end result is more numeric columns will have the same color
gradient (blue-white-red).
* `nmatlist2heatmaps()` displays legend with distinct color steps
for `anno_df` columns with 10 or fewer distinct values, which is
mostly helpful with a small number of integer values.
* `nmatlist2heatmaps()` uses `nrow` for discrete annotation
color layout instead of `ncol` -- jokergoo::ComplexHeatmap 
kindly updated the package to handle both, but nrow is more
definitive in this use case.

# platjam 0.0.25.900

## bug fixes

* `applyXlsxConditionalFormatByColumn()` fixed bug in verbose
output that referred to nonexistent object. Also removed `verbose=TRUE`
from calling function `save_salmon_qc_xlsx()`.

# platjam 0.0.24.900

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

# platjam 0.0.23.900

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

# platjam 0.0.22.900

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

# platjam 0.0.21.900

## changes to existing functions

* `parse_ucsc_gokey()` was updated to recognize bigBed/bb files
can not become overlay tracks.

# platjam 0.0.20.900

## changes to existing functions

* `get_track_defaults()` was modified to include `shortLabel,longLabel`
in the template for multiWig tracks. Apparently a label is displayed
and it uses the internal hub number as a prefix by default.

# platjam 0.0.19.900

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

# platjam 0.0.18.900

## Changes to existing functions

* `get_numeric_transform()` was updated to clean up the help
docs, to accept `"linear"` as a valid name (which does not change
the input data), and to parse the input `transform` values in
order to assign reasonable names where possible.

# platjam 0.0.17.900

## Bug fixes

* Fixed various bugs in `nmatlist2heatmaps()`. No doubt there will
be more, covering the host of assumptions made for the panel_groups,
signal_ceiling, and ylims arguments.

# platjam 0.0.16.900

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

# platjam 0.0.15.900

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

# platjam 0.0.14.900

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

# platjam 0.0.13.900

## Enhancements

* `nmatlist2heatmaps()` altered the logic used to define color ranges
in numeric columns of `anno_df`, opting to use quantile ranges
`c(0.005, 0.995)` to help trim extreme values. Also changed `"purple"`
to `"Purples"` to use the `RColorBrewer` purple color gradient.

# platjam 0.0.12.900

## bug fixes

* `get_numeric_transform()` bug fixed wrong argument name used inside
the function. Should also resolve issues with `nmatlist2heatmaps()`.

# platjam 0.0.11.900

## new functions

* `import_nanostring_rcc()` reads Nanostring RCC files, and produces
a `SummarizedExperiment` object suitable for use in omics-style
analyses.

# platjam 0.0.10.900

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


# platjam 0.0.9.900

## enhancements

* `coverage_matrix2nmat()` now accepts a vector of files as
input, and will process them and return a list of normalizedMatrix
objects, ready to go directly into `nmatlist2heatmaps()`.

# platjam 0.0.8.900

## bug fixes

* `nmatlist2heatmaps()` fixed bug in placement of argument
show_error, now properly placed inside `EnrichedHeatmap::anno_enriched()`.

# platjam 0.0.7.900

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

# platjam 0.0.6.900

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

# platjam 0.0.5.900

## changes

* Version bump while debugging logistics with 
`nmatlist2heatmaps()` argument `transform`.

# platjam 0.0.4.900

## changes

* `nmatlist2heatmaps()` argument `transform` is now applied
to each heatmap, with `transform` values recycled to the
`length(nmatlist)`, to allow control over mathematical
transformation of each numeric matric.

# platjam 0.0.3.900

## changes

* `nmatlist2heatmaps()` argument anno_df is a data.frame with
annotations that can be used to order rows in the data.
* `nmatlist2heatmaps()` updated its method to handle proper
order of user-supplied partition information.
* `nmatlist2heatmaps()` new arguments to handle row labeling,
still testing the combination of partitions and row labels
which appear misaligned.

# platjam 0.0.2.900

## changes

* The package itself allows installation on R-3.0.0 and higher,
resolving an issue with installing on R-3.5.0 when the
package said it required R-3.5.1 (and did not).

# platjam 0.0.2.900

## new function

* `nmatlist2heatmaps()` creates multiple coverage heatmaps
as the result of importing data with `coverage_matrix2nmat()`.

# platjam 0.0.1.900

Initial release.

## new functions

* `coverage_matrix2nmat()` imports genome coverage matrix
file and converts to class `normalizedMatrix` for use in
the amazing `"EnrichedHeatmap"` package.
