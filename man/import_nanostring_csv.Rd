% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/platjam-import-nanostring-csv.R
\name{import_nanostring_csv}
\alias{import_nanostring_csv}
\title{Import Nanostring csv file}
\usage{
import_nanostring_csv(
  csv,
  probe_colname = "Probe_ID",
  probe_anno_file = NULL,
  ...
)
}
\arguments{
\item{csv}{\code{character} path to CSV file, suitable for import
by \code{data.table::fread()}, so most delimiters will be recognized.}

\item{probe_colname}{\code{character} string used to define the probe
identifier in \code{rowData()}, derived from the \code{colnames()} of the
import data. Sometimes the values represent \code{"Probe_ID"}, and
sometimes they may represent \code{"Gene"}, \code{"miRNA"} or another
measurement.}

\item{probe_anno_file}{\code{character} path to one or more files with
additional annotation for each measurement, with identifiers that
match \code{probe_colname} above.}

\item{...}{additional arguments are ignored.}
}
\value{
\code{SummarizedExperiment}
}
\description{
Import Nanostring csv file
}
\seealso{
Other jam import functions: 
\code{\link{coverage_matrix2nmat}()},
\code{\link{deepTools_matrix2nmat}()},
\code{\link{frequency_matrix2nmat}()},
\code{\link{import_lipotype_csv}()},
\code{\link{import_nanostring_rcc}()},
\code{\link{import_proteomics_PD}()},
\code{\link{import_salmon_quant}()}
}
\concept{jam import functions}