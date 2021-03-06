% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterluts.R
\name{col.rep}
\alias{col.rep}
\title{color repetition}
\usage{
col.rep(col, n, alpha = FALSE)
}
\arguments{
\item{col}{color (valid input to \code{\link[grDevices]{col2rgb}}}

\item{n}{number of repetitions (should be a positive integer value)}

\item{alpha}{either logical or actual alpha (in [0..255])}
}
\value{
color matrix with \code{n} columns of color \code{col}
}
\description{
Repeat color \code{col} \code{n} times as columns in a \eqn{3\times n}{3xn}
RGB color matrix. Analog to \code{\link[base]{rep}} which is used internally.
}
\details{
If \code{col} consists of more than one color, those colors are recycl
the usual rules to reach \code{n} colors. In the special case where \code{n}
is smaller than thenumber of colors in \code{col}, only the first \code{n}
colors will be returned.

For \code{n}<1 an empty matrix is returned.
}
\examples{
col.rep('red', 3)
}
\author{
Benno Pütz \email{puetz@psych.mpg.de}
}
