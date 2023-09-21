
#' Validate coverage heatmap params for Rmarkdown
#'
#' Validate coverage heatmap params for Rmarkdown (INCOMPLETE)
#'
#' @family jam coverage heatmap functions
#'
#' @export
validate_heatmap_params <- function
(param,
 ...)
{
   # coverage_file
   coverage_file <- param$coverage_file;
   if ("data.frame" %in% class(coverage_file)) {
      # use as-is
   } else {
      if (!file.exist(coverage_file)) {
         stop(paste0("coverage_file file is not accessible:\n'",
            coverage_file, "'."))
      }
      coverage_file <- data.table::fream(coverage_file,
         data.table=FALSE);
   }
   # validate coverage_file colnames
   # filename
   filename_colname <- head(provigrep(c("^file", "^filename", "^file.name", "file"),
      colnames(coverage_file)), 1);
   if (length(filename_colname) == 0) {
      stop("Input coverage_file did not contain a column name with 'file' or 'filename'.");
   }
   # label
   if (!"label" %in% colnames(coverage_file)) {
      # make unique label based upon the base filename
      coverage_file$label <- gsub("[.][^.]*$",
         "",
         basename(coverage_file[[filename_colname]]));
   }
   # code - used mostly to define contrasts
   if (!"code" %in% colnames(coverage_file)) {
      coverage_file$code <- jamba::colNum2excelName(seq_len(nrow(coverage_file)));
   } else {
      # convert all non-alphanumeric to underscore
      coverage_file$code <- gsub("[^A-Za-z0-9]+",
         "_",
         coverage_file$code);
   }
   # check for files that exist, otherwise validate contrasts
   iexist <- file.exists(coverage_file[[filename_colname]]);
   if (any(!iexist)) {
   }

   # panel_group
   if (!"panel_group" %in% colnames(coverage_file)) {
      panel_group <- NULL;
   } else {
      panel_group <- coverage_file[,"param_group"]
   }
   # panel_group
   if (!"color_ramp" %in% colnames(coverage_file)) {
      color_ramp <- NULL;
   } else {
      panel_group <- coverage_file[,"param_group"]
   }
}
