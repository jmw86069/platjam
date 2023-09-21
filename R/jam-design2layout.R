
#' Convert experiment design into reasonable plot panel layout
#'
#' Convert experiment design into reasonable plot panel layout
#'
#' This function is in development, and is intended for a specific
#' scenario where an experiment design with optional batch/run
#' can be represented in a sensible layout, with appropriate
#' number of panel columns, and number of panel rows.
#'
#' @family jam utility functions
#'
#' @returns `list` with elements:
#' * `nrow`: number of panel rows
#' * `ncol`: number of panel columns
#' * `blank_pos`: integer vector for blank panels, if any
#' * `Group_ctm`: summary info about the panel layout
#' * `isamples`: character vector indicating order to plot
#'
#' @export
design2layout <- function
(SE,
 group_colname=NULL,
 run_colname=NULL,
 ...)
{
   blank_pos <- NULL;
   if (!"data.frame" %in% class(SE)) {
      SEdf <- data.frame(check.names=FALSE,
         SummarizedExperiment::colData(SE));
   } else {
      SEdf <- SE;
   }
   if (length(run_colname) == 0) {
      run_colname <- "temprun";
      SEdf[[run_colname]] <- "temprun";
   }
   ma_ncol <- ceiling(sqrt(nrow(SEdf)));
   ma_nrow <- ceiling(nrow(SEdf) / ma_ncol);
   Group_ctm <- NULL;

   # iterate each column and convert to factor if needed
   all_colnames1 <- gsub("^[-]", "",
      c(group_colname,
         run_colname));
   for (xcol in all_colnames1) {
      if (!is.factor(SEdf[[xcol]])) {
         SEdf[[xcol]] <- factor(SEdf[[xcol]],
            levels=unique(SEdf[[xcol]]));
      }
   }
   group_values <- jamba::pasteByRowOrdered(SEdf[,group_colname, drop=FALSE]);
   run_values <- jamba::pasteByRowOrdered(SEdf[,run_colname, drop=FALSE]);
   SEdf1 <- SEdf[order(group_values, run_values),, drop=FALSE]

   tryCatch({
      Group_ctm <- as.matrix(table(
         SEdf[,c(group_colname, run_colname), drop=FALSE]));
      Group_ctm <- as.matrix(table(
         group_values,
         run_values));
      sub_n <- apply(Group_ctm, 2, function(i){
         max(i, na.rm=TRUE)
      });
      ctm_n <- matrix(cumsum(rep(sub_n, nrow(Group_ctm))),
         ncol=ncol(Group_ctm),
         byrow=TRUE);
      #ctm_n;
      ma_ncol <- sum(sub_n);
      ma_nrow <- nrow(Group_ctm);
      blank_pos <- unlist(lapply(seq_len(nrow(Group_ctm)), function(i){
         jamba::rmNULL(lapply(seq_len(ncol(Group_ctm)), function(j){
            iseq <- seq(from=ctm_n[i,j], by=-1, length.out=sub_n[j]);
            head(iseq, sub_n[j] - Group_ctm[i,j])
         }))
      }));
   }, error=function(e){
      print(e);
   });
   list(nrow=ma_nrow,
      ncol=ma_ncol,
      blank_pos=blank_pos,
      Group_ctm=Group_ctm,
      isamples=rownames(SEdf1))
}
