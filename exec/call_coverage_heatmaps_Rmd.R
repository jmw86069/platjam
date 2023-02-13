## Rscript wrapper to Coverage_Heatmaps_07dec2022.Rmd
#
# This script should probably be called:
# call_coverage_heatmaps_Rmd.Rscript

# pattern to render Rmarkdown with parameters
rmarkdown::render("MyDocument.Rmd", params = list(
   year = 2017,
   region = "Asia",
   printcode = FALSE,
   file = "file2.csv"
))

