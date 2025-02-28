
# geomx_functions.R
#
# Functions to help manipulate and manage GeoMx related data processing


#' Import GeoMx SampleSheet.csv data
#'
#' Import GeoMx SampleSheet.csv data
#'
#' This function can import `"SampleSheet.csv"` and `"GNP_config.ini"`
#' files, which are characterized as follows:
#'
#' * Each subset of data is preceded by a header line: `[header_name]`
#' * Data following this line is comma-delimited, or delimited with
#' `" = "`, both of which are treated as equivalent.
#' * There is a blank line between the subset of data and the next header.
#'
#' Some rules with `return_type="dflist"`:
#'
#' * Delimiters are recognized as `" = "` or `","`, but because the
#' import process calls `data.table::fread()` it will probably also
#' accept tab-delimited data.
#' * When the first row following the header appears to have column names,
#' they are used as-is as column names.
#' * The following criteria cause the first row NOT to be used as
#' column header:
#'
#'    * The first entry does not contain `"Sample_ID"`, and any of:
#'    * Any value following a comma begins with a number, or
#'    * Any value following a comma is `"true"` or `"false"`, or
#'    * Any value following a comma is purely DNA sequence `"[ATGC]+"`, or
#'    * There is no `","` delimiter, or
#'    * There is only one value in the subset of data.
#'
#' Therefore, when the first entry begins with `"Sample_ID"` the
#' first entry is used as column header.
#' * When the first entry is not used as column headers, the heading name
#' itself is used as the first column name, followed by `V` concatenated
#' to the integer column number, for example: `"Sequencing", "V2", "V3", "V4"`
#'
#' ## Todo:
#'
#' * Add `return_type` option to create commands to rename demux output
#' files to the expected GeoMx Sample_ID format.
#'
#' @family jam GeoMx functions
#'
#' @returns `list` or `character` vector, consistent with `return_type`.
#'
#' @param x `character` path to SampleSheet.csv formatted GeoMx file.
#' @param return_type `character` string to define the return type:
#' * `"dflist"` - `list` of `data.frame` for each heading (default)
#' * `"list"` - `list` of lines as-is.
#' * `"indices"` - `character` vector of expected indices, in format
#'    `"Index+Index2"`.
#' * `"filenames"` - `data.frame` describing the expected output filenames
#'    after running `demuxbyname.sh` (BBTools), and the expected
#'    GeoMx filename using the Sample_ID, suitable to rename one file
#'    to the other.
#' @param do_revcomp `logical` indicating whether to use reverse complement
#'    for Index2, default=TRUE.
#' @param demux_prefix `character` prefix assigned when running
#'    `demuxbyname.sh` (BBTools). It is usually in the form
#'    `"project_hist1"` where `"hist1"` refers to the Hamming distance
#'    threshold `"hdist"` used with `demuxbyname.sh`.
#'    This prefix therefore is used to match the filename produced by
#'    `demuxbyname.sh`, which is then renamed to the NGS filename.
#' @param demux_sheetnumber `character` used when `return_type="filenames"`
#'    in formulating the correct GeoMx input filename for the NGS pipeline.
#'    This value therefore affects the filename after renaming the demux
#'    file, to become the NGS input filename.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' samplefile <- system.file("data", "SampleSheet.csv", package="platjam")
#' samplelist <- import_geomx_samplesheet(samplefile)
#' lengths(samplelist)
#'
#' @export
import_geomx_samplesheet <- function
(x="SampleSheet.csv",
 return_type=c("dflist",
    "list",
    "indices",
    "filenames"),
 do_revcomp=TRUE,
 demux_sheetnumber="S1",
 demux_prefix="BK-Gq-1_hdist1_",
 ...)
{
   #
   return_type <- match.arg(return_type)
   #
   if (length(x) == 0) {
      stop("x must be the path to a GeoMx file, which was empty.")
   }
   x <- head(x, 1);
   if (!file.exists(x)) {
      stop("x must be the path to GeoMx file, file not found.")
   }
   samplelines <- readLines(x)

   hlines <- grep("^[[]", samplelines)
   hnames <- gsub(",", "", gsub("[[]|]", "", samplelines[hlines]))
   samplelines_hdf <- data.frame(check.names=FALSE,
      row.names=hnames,
      header=hnames,
      a1=hlines + 1,
      a2=c(tail(hlines - 1, -1), length(samplelines)))

   samplelist <- lapply(jamba::nameVector(samplelines_hdf$header), function(i){
      a1 <- samplelines_hdf[i, "a1"];
      a2 <- samplelines_hdf[i, "a2"];
      v <- samplelines[seq(from=a1, to=a2)];
      v <- v[!v %in% c("", NA)];
      if (all(grepl("=", v))) {
         v <- gsub(" = ", ",", v)
      }
      if (!"list" %in% return_type) {
         if (!grepl("^Sample_ID($|,)", v[1]) && (
               length(v) == 1 ||
               grepl(",([-.0-9]|(true|false)($|,))", v[1]) ||
               grepl(",[ATGC]+($|,)", v[1]) ||
               !grepl(",", v[1]))) {
            vdf <- data.table::fread(
               header=FALSE,
               text=v,
               data.table=FALSE)
            colnames(vdf) <- paste0("V", seq_len(ncol(vdf)));
            colnames(vdf)[1] <- i;
            vdf;
         } else {
            vdf <- data.table::fread(
               text=v,
               data.table=FALSE)
            vdf;
         }
      } else {
         v;
      }
   })

   if ("BCLConvert_Data" %in% names(samplelist) &&
         !"list" %in% return_type) {
      index_colname <- head(jamba::vigrep("index$",
         colnames(samplelist$BCLConvert_Data)), 1)
      index2_colname <- head(jamba::vigrep("index2$",
         colnames(samplelist$BCLConvert_Data)), 1)
      if (length(index_colname) == 0 || length(index2_colname) == 0) {
         stop("BCLConvert_Data must contain colnames 'Index' and 'Index2'")
      }
      use_index2 <- index2_colname;
      if (TRUE %in% do_revcomp) {
         index2_revcomp_colname <- paste0(index2_colname, "_revcomp");
         samplelist$BCLConvert_Data[[index2_revcomp_colname]] <- revcomp(
            samplelist$BCLConvert_Data[[index2_colname]])
         use_index2 <- index2_revcomp_colname;
      }
      if ("indices" %in% return_type) {
         indices <- paste0(
            samplelist$BCLConvert_Data[[index_colname]],
            "+",
            samplelist$BCLConvert_Data[[use_index2]])
         return(indices)
      } else if ("filenames" %in% return_type) {
         #
         bcldf <- samplelist$BCLConvert_Data;
         demux_base <- paste0(demux_prefix,
            bcldf[[index_colname]],
            "+",
            bcldf[[use_index2]])
         demux_files <- paste0(rep(demux_base, each=2),
            "_", c(1, 2),
            ".fq.gz")
         demux_samplefile <- paste0(
            rep(bcldf[["Sample_ID"]], each=2),
            "_", demux_sheetnumber,
            "_L001_R", c(1, 2),
            "_001.fastq.gz")
         bcldf2 <- data.frame(check.names=FALSE,
            bcldf,
            file_from=matrix(demux_files, ncol=2, byrow=TRUE),
            file_to=matrix(demux_samplefile, ncol=2, byrow=TRUE))
         return(bcldf2)
      }
   } else if (any(c("filenames", "indices") %in% return_type)) {
      stop(paste0(
         "BCLConvert_Data heading was not found, ",
         "required for return_type 'indices' or 'filenames'"))
   }
   return(samplelist)
}


#' Simple DNA reverse-complement function
#'
#' @family jam GeoMx functions
#'
#' @returns `character` string with reverse-complemented sequences.
#'
#' @param x `character` string with DNA or RNA sequence, expected to
#'    contain only A,T,G,C,N, and possibly U.
#' @param molecule `character` string with `"DNA"` (default) or `"RNA"`.
#'    Input recognizes either T or U, but will only produce T for DNA,
#'    and U for RNA.
#'
#' @examples
#' dna_in <- "GTTACTCT";
#' revcomp(dna_in)
#'
#' @export
revcomp <- function
(x,
 molecule=c("DNA",
    "RNA"),
 ...)
{
   molecule <- match.arg(molecule);
   # operate only on unique elements
   xu <- unique(x);
   xs <- strsplit(xu, "")
   # convert to complementary base
   from_to <- c(A="T", U="A", T="A", G="C", C="G", N="N")
   if ("RNA" %in% molecule) {
      from_to <- c(A="U", U="A", T="A", G="C", C="G", N="N")
   }
   xu_revcomp <- jamba::cPaste(sep="", lapply(xs, function(ix){
      rev(unname(from_to[toupper(ix)]))
   }))
   xu_revcomp[match(x, xu)]
}

