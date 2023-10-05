
#' RMarkdown tab iterator
#'
#' RMarkdown tab iterator
#'
#' This function is intended to help automate the process of generating
#' tabs in RMarkdown output, particularly when there may be multiple
#' layers of tabs. The driving example plots a figure at the lowest level
#' of tabs, using the options defined by the tabs collectively.
#'
#' While tabs can be created using `for` loops or `lapply()` methods,
#' the RMarkdown code quickly becomes complicated and difficult to maintain.
#'
#' Each layer of tabs is iterated, and a preparatory function can be called,
#' for example to update data to be used downstream.
#' Any layer is free to create a figure, table, or text output, however
#' it is usually recommended to create output at the lowest layer.
#'
#' Output can be defined in the final `fn_list` layer, however it is cleaner
#' to use the argument `base_fn` to generate final output such as figure
#' or table. The functions in `fn_list` should be limited to those required
#' for updates during the process.
#'
#' @returns `NULL` invisibly, this function is called for the by-product
#'    of printing RMarkdown-compatible output.
#'
#' @family jam utility functions
#'
#' @param tab_list `list` of `character` or `list` objects.
#'   * each `list` entry represents a layer of tabs
#'   * `names(tab_list)` defines the variable name to which each element
#'   in the `list` is assigned.
#'   For example `tab_list=list(use_gene=c("ACTB", "GAPDH")` will create
#'   a layer of tabs, with tabs `"ACTB"`, and `"GAPDH"`. For the `"ACTB"`
#'   tab it will assign `use_gene <- "ACTB"` for use in the tab.
#'   * the names of the `list` should match `names(prep_list)`
#'   * each `character` vector contains the tabs to display
#' @param fn_list `list` of `functions` to apply to corresponding
#'    layers of tabs.
#'    * `names(prep_list)` must match `names(tab_list)`, and provides
#'    the function to be applied at the specified layer of tabs.
#'    This function should create the tab content, and should not print
#'    the tab header itself.
#' @param final_fn `function` called at the bottom layer of tabs, after
#'    the corresponding function from `fn_list` is called, if it exists.
#'    This function is optional, and is intended to be a convenient place
#'    for the core function in the stack, without having to include this
#'    function in the final layer of `fn_list`.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' tab_labels <- list(
#'    assay_name=c(
#'       raw="raw data",
#'       jammanorm_raw="log ratio normalized",
#'       quantile_raw="quantile normalized"),
#'    centerby_name=c(
#'       global="global-centered",
#'       `1st_control`="centered vs first control")
#'    )
#' base_fn <- function(...){
#'    plot(x=c(0, 3), y=c(0, 3), asp=1, pch=".")
#'    text(x=c(1, 1.5, 2), y=c(1, 1.5, 2),
#'       c(paste0("assay_name:\n", assay_name),
#'       paste0("formatting:\n", formatting),
#'       paste0("centerby_name:\n", centerby_name)))
#' }
#' tab_list <- list(
#'    assay_name=c("raw", "jammanorm_raw", "quantile_raw"),
#'    formatting=c("tab-delimited", "RData"),
#'    centerby_name=c("global", "1st_control"))
#' rmd_tab_iterator(
#'    tab_list=tab_list,
#'    fn_list=NULL,
#'    base_fn=base_fn,
#'    tab_labels=tab_labels,
#'    heading_level=3,
#'    verbose=FALSE)
#'
#' @export
rmd_tab_iterator <- function
(tab_list,
 fn_list=NULL,
 base_fn=NULL,
 tab_labels=NULL,
 heading_level=2,
 verbose=FALSE,
 htmlOut=TRUE,
 envir=NULL,
 ...)
{
   #
   printRmd <- function(..., htmlOut1=htmlOut, comment=FALSE, timeStamp=FALSE){
      jamba::printDebug(...,
         htmlOut=htmlOut1,
         timeStamp=timeStamp,
         comment=comment)
   }
   # create a specific environment to hold iterated variables
   if (!"environment" %in% class(envir)) {
      if (verbose) {
         printRmd("rmd_tab_iterator(): ",
            "created new environment ", "envir")
      }
      envir <- new.env();
   }
   # attach(envir)
   # suppressWarnings(attach(envir))

   tab_values <- tab_list[[1]];
   tab_name <- head(names(tab_list), 1);
   tab_label_set <- NULL;
   if (tab_name %in% names(tab_labels)) {
      tab_label_set <- tab_labels[[tab_name]];
   }

   heading_string <- paste0(rep("#", heading_level), collapse="");
   for (tab_value in tab_values) {
      if (tab_value %in% names(tab_label_set)) {
         tab_label <- jamba::cPaste(sep=" ", tab_label_set[[tab_value]]);
      } else {
         tab_label <- jamba::cPaste(sep=" ", tab_value);
      }
      tab_suffix <- "{.tabset}";
      if (length(tab_list) == 1) {
         tab_suffix <- "";
      }
      # print the tab header
      # TODO: suppress the header when no content is present,
      # which probably requires assembling each sub-component to determine
      # whether there is any output... Not going to happen for now.
      cat(paste0("\n\n",
         heading_string, " ",
         tab_label, " ",
         tab_suffix,
         "\n\n"))

      # assign value to the environment
      assign(x=tab_name,
         value=tab_value,
         envir=envir)

      if (verbose) {
         printRmd("rmd_tab_iterator(): ",
            tab_name, " <- ",
            c("c(",
               jamba::cPaste(paste0("'", tab_value, "'"), sep=", "),
               ")"),
            sep="")
         # "tab_name: ", tab_name,
         # ", tab_value: ", tab_value);
      }

      # perform tab function
      if (tab_name %in% names(fn_list)) {
         tab_fn <- fn_list[[tab_name]];
         environment(tab_fn) <- envir;
         # call this function
         if (verbose) {
            printRmd("rmd_tab_iterator(): ",
               "calling the tab function");
         }
         tab_fn(...)
      }

      # define the next iterator
      new_tab_list <- tail(tab_list, -1);
      if (length(new_tab_list) > 0) {
         new_fn_list <- fn_list[!names(fn_list) %in% tab_name]
         new_tab_labels <- tab_labels[!names(tab_labels) %in% tab_name]
         new_heading_level <- heading_level + 1;
         # define function so it inherits this enclosing environment
         next_iterator <- function
         (tab_list=new_tab_list,
            fn_list=new_fn_list,
            base_fn1=base_fn,
            tab_labels=new_tab_labels,
            heading_level=new_heading_level,
            htmlOut1=htmlOut,
            verbose1=verbose,
            ...)
         {
            rmd_tab_iterator(tab_list=tab_list,
               fn_list=fn_list,
               base_fn=base_fn1,
               tab_levels=tab_levels,
               heading_level=heading_level,
               htmlOut=htmlOut1,
               verbose=verbose1,
               envir=envir,
               ...)
         }
         if (verbose) {
            printRmd("rmd_tab_iterator(): ",
               "calling the next layer of tabs");
         }
         next_iterator(...)
      } else if ("function" %in% class(base_fn)) {
         # if present, call the final function
         environment(base_fn) <- envir;
         # call this function
         if (verbose) {
            printRmd("rmd_tab_iterator(): ",
               "calling the ", "base_fn", " function");
         }
         base_fn(...)
      }
   }
   return(invisible(NULL))
}
