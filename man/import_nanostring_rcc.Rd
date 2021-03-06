% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/platjam-nanostring.R
\name{import_nanostring_rcc}
\alias{import_nanostring_rcc}
\title{Import Nanostring RCC files}
\usage{
import_nanostring_rcc(
  rcc_files = NULL,
  rcc_path = ".",
  rcc_pattern = NULL,
  exclude = NULL,
  include = NULL,
  nprobes = -1,
  control_greps = c(POS = "^POS_", NEG = "^NEG_"),
  return_type = c("SummarizedExperiment", "NanoString"),
  debug = FALSE,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{rcc_files}{character vector of file paths to RCC files,
or if \code{NULL} then the \code{rcc_path} and \code{rcc_pattern} arguments
are used to find RCC files.}

\item{rcc_path}{character vector of file directories to look
for RCC files, the files are matched to pattern \code{rcc_pattern}.
Note the \code{rcc_path} is only used when \code{rcc_files} is not supplied.}

\item{rcc_pattern}{character string containing regular expression
pattern, used to match filenames which should be considered RCC
files.}

\item{exclude}{optional character vector of filenames to exclude from
import. The filenames may be either a full file path, or
the basename without the file path.}

\item{include}{optional character vector of filenames to include
during data import. Note that when \code{include} is defined, \emph{only}
these files are included in the file import, all other files are
ignored.}

\item{debug}{logical indicating whether to send intermediate data
before full processing, useful for debugging file format errors.}

\item{verbose}{logical indicating whether to print verbose output.}

\item{...}{additional arguments are ignored.}
}
\value{
\code{SummarizedExperiment} or \code{NanoString} object, based
upon argument \code{return_type}.
}
\description{
Import Nanostring RCC files
}
\details{
This function reads one or more Nanostring RCC files,
and produces a \code{SummarizedExperiment} object.

This function was based upon the function
\code{NanoStringNorm::read.markup.RCC()}, which in our experience
was problematic because it assumed all files contained identical
sample annotation content.

For examples how to curate sample names into proper \code{data.frame}
of experimental factors, see \code{splicejam::curateVtoDF()} or
\code{splicejam::curateDFtoDF()}.
}
\seealso{
Other jam import functions: 
\code{\link{coverage_matrix2nmat}()},
\code{\link{deepTools_matrix2nmat}()},
\code{\link{frequency_matrix2nmat}()}
}
\concept{jam import functions}
