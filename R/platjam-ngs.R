
#' Get Salmon metadata and aux info into a data.frame
#'
#' Get Salmon metadata and aux info into a data.frame
#'
#' This function takes a file path to one or more Salmon
#' output files, uses that path to locate the full set of
#' available files, loads data from each of the discovered
#' files, and returns the results in a `data.frame` format.
#'
#' This function uses `rprojroot::find_root()` to find the
#' root directory, defined as the directory that contains
#' the file `"cmd_info.json"`. The path to `"meta_info.json"`
#' is constructed relative to that location.
#'
#' Recognized files:
#'
#' * `meta_info.json` - typically in a subdirectory `aux_info/meta_info.json`
#' * `cmd_info.json` - typically in the same directory as the `aux_info`
#' directory.
#'
#' If a relative path to `"cmd_info.json"` cannot be determined, this
#' function returns `NULL`.
#'
#' When the input `metafile` includes multiple files, only
#' the unique Salmon root directories are returned.
#'
#' This function uses `jsonlite` to read each JSON file, which
#' is converted to a `data.frame`. Any JSON fields that contain
#' multiple values are comma-delimited using `jamba::cPaste()`
#' in order to fit on one row in the `data.frame`.
#'
#' @family jam nextgen sequence functions
#'
#' @return `data.frame` whose number of rows is equal to the number
#'    of unique Salmon root directories in the input `metafile`.
#'    For any input `metafile` not found, the output is `NULL`.
#'
#' @param metafile character vector of one or more files, usually the
#'    full file path to the `meta_info.json` file after running Salmon
#'    quant. The path `metafile` should be the path to any output file
#'    from one Salmon quant analysis.
#' @param exclude_hashes logical indicating whether to drop columns that
#'    contain file hashes.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' cmdinfopath <- system.file("data", "salmonOut", "cmd_info.json", package="platjam");
#' if (nchar(cmdinfopath) > 0) {
#'    get_salmon_meta(cmdinfopath);
#' }
#'
#' @export
get_salmon_meta <- function
(metafile,
 exclude_hashes=TRUE,
 ...)
{
   ##
   if (!suppressPackageStartupMessages(require(rprojroot))) {
      stop("The rprojroot package is required for get_salmon_meta().");
   }
   if (!suppressPackageStartupMessages(require(jsonlite))) {
      stop("The jsonlite package is required for get_salmon_meta().");
   }
   if (length(metafile) > 1) {
      metapaths <- unique(jamba::rmNA(get_salmon_root(metafile)));
      metajsons <- lapply(seq_along(metapaths), function(i){
         get_salmon_meta(metapaths[i],
            exclude_hashes=exclude_hashes,
            ...);
      });
      metajson <- jamba::mergeAllXY(metajsons);
      return(metajson);
   }

   ## Try to use rprojroot to find the Salmon root directory
   metapath <- get_salmon_root(metafile);
   if (is.na(metapath)) {
      return(NULL);
   }
   metafile <- file.path(metapath, "aux_info",
      "meta_info.json");
   cmdinfofile <- file.path(metapath,
      "cmd_info.json");
   libfile <- file.path(metapath,
      "lib_format_counts.json");

   json1 <- jsonlite::read_json(metafile);
   jsoninfo <- jamba::sdim(json1);
   if (any(jsoninfo$class %in% "list" & jsoninfo$rows > 0)) {
      jsonlist <- which(jsoninfo$class %in% "list" & jsoninfo$rows > 0);
      json1[jsonlist] <- jamba::cPaste(json1[jsonlist]);
   }
   if (length(names(metafile)) > 0) {
      metajson <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         name=names(metafile),
         as.data.frame(json1[jsoninfo$rows > 0]));
   } else {
      metajson <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         as.data.frame(json1[jsoninfo$rows > 0]));
   }
   if (exclude_hashes) {
      metajson <- metajson[,jamba::unvigrep("hash", colnames(metajson)),drop=FALSE];
   }

   ## cmd_info
   if (file.exists(cmdinfofile)) {
      json1 <- jsonlite::read_json(cmdinfofile);
      jsoninfo <- jamba::sdim(json1);
      if (any(jsoninfo$class %in% "list" & jsoninfo$rows > 0)) {
         jsonlist <- which(jsoninfo$class %in% "list" & jsoninfo$rows > 0);
         json1[jsonlist] <- jamba::cPaste(json1[jsonlist]);
      }
      if (any(jsoninfo$class %in% "list" & jsoninfo$rows == 0)) {
         jsonlist0 <- which(jsoninfo$class %in% "list" & jsoninfo$rows == 0);
         json1[jsonlist0] <- sapply(names(json1[jsonlist0]), function(i){i});
      }
      cmdjson <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         as.data.frame(json1));
      metajson[,colnames(cmdjson)] <- cmdjson;
   }

   ## lib_format_counts
   if (file.exists(libfile)) {
      json1 <- jsonlite::read_json(libfile);
      jsoninfo <- jamba::sdim(json1);
      if (any(jsoninfo$class %in% "list" & jsoninfo$rows > 0)) {
         jsonlist <- which(jsoninfo$class %in% "list" & jsoninfo$rows > 0);
         json1[jsonlist] <- jamba::cPaste(json1[jsonlist]);
      }
      if (any(jsoninfo$class %in% "list" & jsoninfo$rows == 0)) {
         jsonlist0 <- which(jsoninfo$class %in% "list" & jsoninfo$rows == 0);
         json1[jsonlist0] <- sapply(names(json1[jsonlist0]), function(i){i});
      }
      libjson <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         as.data.frame(json1[jsoninfo$rows > 0]));
      metajson[,colnames(libjson)] <- libjson;
   }

   return(metajson);
}

#' Get Salmon root directory
#'
#' Get Salmon root directory
#'
#' This function uses the `rprojroot` package to find the
#' root output directory of a Salmon quant analysis, looking
#' for the file `"cmd_info.json"`. If this file is not found,
#' this function returns `NULL`.
#'
#' @return file path to the Salmon root directory, or `NA` if
#'    the Salmon root directory could not be found.
#'
#' @param file character path to a Salmon output file, or a vector
#'    of files.
#'
#' @family jam nextgen sequence functions
#'
#' @examples
#' cmdinfopath <- system.file("data", "salmonOut", "cmd_info.json", package="platjam");
#' if (nchar(cmdinfopath) > 0) {
#'    get_salmon_root(cmdinfopath);
#' }
#'
#' @export
get_salmon_root <- function
(file)
{
   ## Try to use rprojroot to find the Salmon root directory
   if (!suppressPackageStartupMessages(require(rprojroot))) {
      stop("The rprojroot package is required for get_salmon_meta().");
   }
   if (length(file) > 1) {
      salmonroots <- sapply(file, get_salmon_root);
      return(salmonroots);
   }
   salmonroot <- tryCatch({
      rprojroot::find_root(path=file,
         "cmd_info.json");
   }, error=function(e){
      NA;
   });
   return(salmonroot);
}
