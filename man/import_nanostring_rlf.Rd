% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/platjam-nanostring-rlf.R
\name{import_nanostring_rlf}
\alias{import_nanostring_rlf}
\title{Import Nanostring RLF file for a codeset}
\usage{
import_nanostring_rlf(
  rlf,
  plot_type = c("none", "table", "list"),
  color_sub = c(Endogenous = "darkorange", Positive = "greenyellow", Negative =
    "firebrick3", Housekeeping = "gold", Binding = "steelblue4", Purification =
    "steelblue2", Reserved = "mediumpurple1", ClassDate = "blue", ClassName_count =
    "navy"),
  desat = c(0.2, 0.4),
  ...
)
}
\arguments{
\item{rlf}{\code{character} path to a Nanostring RLF file}

\item{plot_type}{\code{character} string passed to \code{design2colors()} to assign
colors to RLF data.
Default \code{plot_type="none"} creates no plot.}

\item{color_sub, desat}{arguments passed to \code{design2colors()}, where
\code{color_sub} is a name-color pairing to assign specific colors to
values, or colnames.}

\item{...}{additional arguments are passed to \code{design2colors()}.}
}
\value{
\code{SummarizedExperiment} object
}
\description{
Import Nanostring RLF file for a codeset
}
\details{
This function partially parses and imports data encoded in RLF
file format by Nanostring for a particular Nanostring codeset.
This file often includes information that links each probe to
a "Class", which is then associated with other helpful information:
\itemize{
\item "ClassKey" - usually an integer unique identifier
\item "ClassName" - character string label for each ClassKey
\item "ClassActive" - integer value indicating the activity of probe for
downstream stats.
\item "ClassDate" - the date associated with each class, usually identical
for all classes in one codeset.
\item "ClassSource" - optional field describing the origin of each class
\item "ClassPreparer" - optional tool used to prepare classes in the RLF file
}

The most useful fields for downstream analysis are:
\itemize{
\item "ClassName" - which contains:
\itemize{
\item "Positive", "Negative" for the positive and negative control probes
\item "Endogenous" for the intended biological measurements, usually in
the form of genes, transcripts, or miRNAs.
\item "Housekeeping" for the intended low-variance conrol measurements,
usually in the form of genes, transcripts, miRNAs that are known or
expected to have stable expression across the experiment conditions.
\item "Binding", "Purification", "Reserved" describe technical features
beyond the scope of most analysis work.
}
\item "ClassActive" - which contains:
\itemize{
\item \code{2} - a value that links \code{"Endogenous"} and \code{"Housekeeping"} above,
though this value is not associated with another label described in
the RLF file.
\item \code{1} - a value that links \code{"Binding"}, \code{"Purification"}, \code{"Reserved"}
\item \code{3} - a value that links \code{"Positive"} and \code{"Negative"} as control
spike-in probes on the codeset
}
}
}
\seealso{
Other jam import functions: 
\code{\link{coverage_matrix2nmat}()},
\code{\link{deepTools_matrix2nmat}()},
\code{\link{frequency_matrix2nmat}()},
\code{\link{import_lipotype_csv}()},
\code{\link{import_metabolomics_niehs}()},
\code{\link{import_nanostring_csv}()},
\code{\link{import_nanostring_rcc}()},
\code{\link{import_proteomics_PD}()},
\code{\link{import_proteomics_mascot}()},
\code{\link{import_salmon_quant}()},
\code{\link{process_metab_compounds_file}()}
}
\concept{jam import functions}
