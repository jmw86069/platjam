#' Import deepTools coverage matrix to normalizedMatrix
#'
#' Import deepTools coverage matrix to normalizedMatrix
#'
#' This function is under active development.
#'
#' @family jam import functions
#'
#' @export
deepTools_matrix2nmat <- function
(x=NULL,
 filename=NULL,
 signal_name=NULL,
 target_name="target",
 background=0,
 smooth=FALSE,
 target_is_single_point=FALSE,
 signal_is_categorical=FALSE,
 mat_grep="[-0-9]+:[-0-9]+",
 upstream_grep="^[-]",
 downstream_grep="^[^-]",
 target_grep="^0$",
 verbose=FALSE,
 ...)
{
   x1 <- readLines(filename, n=1);
   xyaml <- yaml::read_yaml(text=gsub("^@", "", x1));
   samples <- xyaml$sample_labels;
   samples;
   binsize <- xyaml$`bin size`;
   binsize;
   upstream <- xyaml$upstream;
   upstream;
   group_boundaries <- xyaml$group_boundaries;
   sample_boundaries <- xyaml$sample_boundaries;
   x <- data.table::fread(text=readLines(filename),
      skip=1,
      header=FALSE,
      sep="\t");
   rownames(x) <- x[[4]];

   mats <- as.matrix(x[,-1:-6,drop=FALSE]);
   mat_n <- ceiling(ncol(mats) / length(samples));

   ## required attributes
   starts <- seq(from=-upstream, by=binsize, length.out=mat_n);
   upstream_index <- which(starts < 0);
   target_index <- which(starts == 0);
   downstream_index <- which(starts > 0);
   target_name <- xyaml$group_labels;

   mat_split <- rep(rep(samples, each=mat_n), length.out=ncol(mats));
   colnames_l <- split(colnames(mats), mat_split);
   mat_l <- lapply(seq_along(colnames_l), function(k){
      i <- colnames_l[[k]];
      im <- mats[,i];
      colnames(im) <- paste0("u", seq_len(ncol(im)));
      rownames(im) <- x[[4]];
      attr(im, "signal_name") <- samples[k];
      attr(im, "target_name") <- target_name;
      attr(im, "upstream_index") <- upstream_index;
      attr(im, "target_index") <- target_index;
      attr(im, "downstream_index") <- downstream_index;
      im;
   });

}
