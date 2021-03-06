% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/platjam-coverhm.R
\name{make_coverage_heatmaps}
\alias{make_coverage_heatmaps}
\title{Make coverage heatmaps}
\usage{
make_coverage_heatmaps(config_df, anno_df = NULL, ...)
}
\description{
(IN DEV) Make coverage heatmaps using a simple set of config options
}
\details{
This function is still in development and is not yet active.

This function is intended as a wrapper function which calls
\code{nmatlist2heatmaps()} in a more organized way.

The \code{config_df} is a \code{data.frame} with one row per coverage heatmap
to produce, and all options associated with that heatmap
are stored on the same row of the \code{data.frame}.
\itemize{
\item file: \code{character} file path to the coverage matrix file
\item name: \code{character} string used as a name, used in difference calculations
\item label: \code{character} string used as a label above each heatmap
\item display: \code{logical} indicating whether to display each heatmap
\item group: \code{character} string indicating a panel group, used to share
certain settings across groups of heatmap panels. When this value
is not defined, each heatmap is defined as its own \code{group}.
\item color: \code{character} string that refers to a specific color gradient,
or a comma-delimited series of colors to use as a color gradient.
See below for details.
\item ceiling: \code{numeric} used to define the maximum numeric value
applied to the color gradient
\item ylim: \code{numeric} used to define a fixed y-axis range used for the
profile plot atop each heatmap. When this value is absent or \code{NA}
the maximum y-axis value for the \code{group} is used.
\item control: \code{character} string that defines an optional control \code{name}
used to subtract coverage from this row.
}
\subsection{Colors}{

Colors can be define one of a few ways:
\itemize{
\item name of a color gradient
\item name of a single color
\item comma-delimited colors, either as color names, or as hex colors
in the format \verb{#FFAABB}.
}

Color gradients from RColorBrewer are recognized. You can view
\itemize{
\item linear color gradients: \code{RColorBrewer::display.brewer.all(type="seq")}
\item divergent color gradients: \code{RColorBrewer::display.brewer.all(type="div")}
}

Color gradients are also recognized from viridis, however these are
only linear: \code{viridis::inferno()}, \code{viridis::magma()}, \code{viridis::plasma()}
\code{viridis::cividis()},\code{viridis::viridis()}.
}
}
\examples{
# display RColorBrewer linear gradients
RColorBrewer::display.brewer.all(type="seq")

# display RColorBrewer divergent gradients
RColorBrewer::display.brewer.all(type="div")

# viridis linear gradients
jamba::showColors(list(
   inferno=viridis::inferno(11),
   magma=viridis::magma(11),
   plasma=viridis::plasma(11),
   cividis=viridis::cividis(11),
   viridis=viridis::viridis(11)))

# jam_linear
showColors(jam_linear)

# jam_divergent
showColors(jam_divergent)

}
