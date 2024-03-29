% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jam-numeric-to-colorfn.R
\name{assign_numeric_colors}
\alias{assign_numeric_colors}
\title{Assign color function to a numeric vector}
\usage{
assign_numeric_colors(
  x,
  color,
  color_max = NULL,
  ratio_for_zero_baseline = 3,
  lower_range_expansion = 0.2,
  base_color = "#fff6f4",
  restrict_pretty_range = TRUE,
  lens = 1,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{x}{\code{numeric} or \code{integer} vector}

\item{color}{\code{character} color}

\item{color_max}{\code{numeric} optional fixed value to define the
\code{numeric} value associated with the maximum color gradient color.}

\item{ratio_for_zero_baseline}{\code{numeric} indicating the ratio of
max to min numeric range, above which the baseline value
should include zero. This argument is intended when values
include \code{c(3, 100)} and the numeric range is more reasonable
to include zero; compared to values between \code{c(185, 192)} where
it makes more sense to apply the color gradient only near
those values, and not to use \code{c(0, 192)}.}

\item{lower_range_expansion}{\code{numeric} used when the baseline is not
zero, see \code{ratio_for_zero_baseline}, to adjust the baseline below
the lowest value by this fraction of the span of values (max - min)
below the lowest value. The intention is to prevent assigning
white \code{base_color} to this lowest value, instead it should use
a very pale variant of \code{color}.}

\item{base_color}{\code{character} string used as the baseline color when
defining a color gradient to \code{color}. This value is passed to
\code{jamba::getColorRamp(..., defaultBaseColor)}.}

\item{restrict_pretty_range}{\code{logical} indicating whether the values
returned by \code{pretty()} should be filtered by \code{range(x)}, excluding
any numeric breaks outside that range. The practical effect is that
it reduces labels displayed by \code{jamba::showColors()} or by color
legends drawn by \code{ComplexHeatmap::Heatmap()}.}

\item{lens}{\code{numeric} value that defines the intensity of color gradient,
with higher values increasing the rate of increase in color intensity,
and negative values decreasing this rate. This argument is passed
to \code{jamba::getColorRamp()}.}

\item{...}{additional arguments are ignored.}
}
\value{
\code{function} as defined by \code{circlize::colorRamp2()} which takes
a \code{numeric} vector as input, and returns \code{character} vector of
colors as output. The attributes described below are used to
show a suitable summary of colors in \code{jamba::showColors()};
and are used by \code{jamses::heatmap_se()} to define a
usable color legend with reasonable number of labels.
The \code{function} also contains two important attributes:
\itemize{
\item \code{"breaks"}: \code{numeric} values at break positions used when the
color function was defined. Note the values are taken from the
color function environment, so modifying breaks directly will
not affect the color function output.
\item \code{"colors"}: a matrix of colors in rgb format, with columns
red, blue, green, and \code{numeric} values ranging from 0 to 1.
The number of rows equals the length of color breaks.
}
}
\description{
Assign color function to a numeric vector
}
\details{
This function is called internally by \code{design2colors()}. It takes
a \code{numeric} vector as input, and applies a color gradient
defined by \code{color} across the range of \code{numeric} values.

When there are negative values, a divergent color scale is generated:
\itemize{
\item When \code{color} is a single R color, it is used as the positive color
gradient. A complementary color is chosen as the negative color
using \code{quick_complement_color()}. The middle color is \code{base_color}.
\item When \code{color} is a named color gradient, it is assumed to be a
divergent gradient. The middle color of the gradient is assigned to zero,
and the color gradient is symmetric above and below zero.
}

The \code{numeric} range is either defined by the data \code{x}, or when
\code{color_max} is supplied, it is used to define the \code{numeric}
value at which (or higher) the maximum color is assigned.
\itemize{
\item When \code{color_max} is not assigned, it is defined as \code{max(abs(x))}.
\item If there are negative values, the color gradient will be
defined from \code{-color_max} to \code{color_max}, with the middle color
assigned to zero, as described above.
\item If there are only positive values, the baseline will be defined
using \code{ratio_for_zero_baseline}, and the maximum color will be defined
at \code{color_max} (and higher).
}

The \code{ratio_for_zero_baseline} is used only when there are
no negative values, since it determines whether the color gradient
should be applied starting at zero, or based upon \code{min(x)}.
When the lowest value is less than \code{ratio_for_zero_baseline * color_max}
the color gradient is applied starting from zero;
otherwise the baseline uses \code{min(x)}. In the latter case, non-zero
baseline, the lowest value is extended by multiplying
\code{lower_range_expansion} by the span of the numeric range, in order
to prevent the lowest non-zero value from being a blank color
defined by \code{base_color}.
}
\seealso{
Other jam color functions: 
\code{\link{design2colors}()},
\code{\link{df_to_numcolors}()},
\code{\link{mean_hue}()},
\code{\link{print_color_list}()},
\code{\link{quick_complement_color}()}
}
\concept{jam color functions}
