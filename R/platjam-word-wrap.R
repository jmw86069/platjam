
#' Apply word-wrap to text, recognizing whitespace and some punctuation
#'
#' @return `character` vector by default, however when `sep=NULL`
#'    it returns a `list` of `character` vectors.
#'
#' @family jam utility
#'
#' @param x `character` string or vector.
#' @param width `numeric` desired maximum character width
#' @param sep `character` string used as a separator between wrapped
#'    subsets, default `"\n"` wraps to new lines.
#'    When `sep=NULL` the wrapped subsets are not combined, instead
#'    are returned as a `list` of `character` vectors.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' x <- c("logFC dH1A_VEH_ATAC-UL3_VEH_ATAC",
#'    "logFC dH1A_VEH_ATAC_UL3_VEH_ATAC",
#'    "one-two-three-four-five-six-seven")
#' cat(apply_word_wrap(x), sep="\n--\n")
#'
#' apply_word_wrap(x, sep=NULL)
#'
#' @export
apply_word_wrap <- function
(x,
 width=15,
 sep="\n",
 ...)
{
   #
   x_list <- lapply(x, function(i){
      # i2 <- gsub("([-,]+)($|[^ ])", "\\1 \\2", i);
      ivals <- stringi::stri_wrap(i,
         width=width,
         whitespace_only=FALSE);
      if (length(sep) == 1) {
         ivals <- jamba::cPaste(ivals, sep=sep)
      }
      ivals
   })
   x_list
   if (length(sep) == 1) {
      return(unlist(x_list))
   }
   x_list
}
