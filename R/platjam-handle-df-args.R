
#' Handle function arguments passed in the form of a data.frame
#'
#' Handle function arguments passed in the form of a data.frame
#'
#' This function is an internal utility function intended to allow
#' funtion arguments to be passed as a single `data.frame`, and
#' each of the `colnames(df)` represent formal function arguments.
#'
#' One tricky aspect is how to handle arguments that only expect
#' one value. It could be validated here, or that could be
#' responsibility of the parent function.
#'
#' @param df `data.frame` with `colnames()` that should be applied to
#'    function arguments defined by `parent_fn`; or a `character`
#'    filename that contains delimited data that can be coerced to
#'    a `data.frame`, for example a tab-delimited or comma-delimited
#'    text file.
#' @param parent_fn `function` used to obtain the `formals()` to define
#'    function arguments.
#' @param ... additional arguments are ignored.
#'
#' @family jam utility functions
#'
#' @returns `environment` containing variables parsed from `df`
#'
#' @examples
#' parent_fn <- nmatlist2heatmaps;
#'
#' @export
handle_df_args <- function
(df,
 exclude=NULL,
 parent_fn=NULL,
 fuzzy_match=TRUE,
 fn_env=NULL,
 ...)
{
   # validate df
   if (length(df) == 0) {
      return(fn_env)
   }
   if (!"data.frame" %in% class(df)) {
      if ("character" %in% df && length(df) == 1) {
         if (grepl("[.]xlsx", ignore.case=TRUE, df)) {
            df <- jamba::readOpenxlsx(df,
               sheet=1,
               ...)[[1]];
         } else {
            df <- data.table::fread(df,
               data.table=FALSE,
               ...);
         }
      } else {
         stop("Input df must be a data.frame, or text file, or .xlsx file");
      }
   }
   # validate fn_env
   if ("NULL" %in% class(fn_env)) {
      fn_env <- new.env();
   }
   #
   parent_args <- setdiff(names(formals(parent_fn)),
      c(exclude, "..."));
   df_args <- colnames(df);
   if (fuzzy_match) {
      arg_adj <- function(x){
         gsub("[^a-z0-9]", ".", tolower(x))
      }
   } else {
      arg_adj <- function(x){
         x
      }
   }
   arg_match <- match(arg_adj(parent_args),
      arg_adj(df_args));
   arg_from <- df_args[arg_match[!is.na(arg_match)]];
   arg_to <- parent_args[!is.na(arg_match)];
   arg_df <- data.frame(arg_from, arg_to);
   print(arg_df)
   for (i in seq_along(arg_to)) {
      assign(envir=fn_env,
         x=arg_to[i],
         value=df[[arg_from[i]]]);
      jamba::printDebug("assign x='", arg_to[i], "'='",
         df[[arg_from[i]]]);
   }
   return(fn_env);
}
