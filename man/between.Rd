% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterluts.R
\name{between}
\alias{between}
\title{between test}
\usage{
between(x, low, high = NULL, index = TRUE, named = TRUE, edges = NULL, ...)
}
\arguments{
\item{x}{value(s) to test}

\item{low}{lower bound(s) of test interval(s)}

\item{high}{upper bound(s) of test interval(s)}

\item{index}{when set return indices, otherwise match matrix (see description)}

\item{named}{whether to put names on the return value}

\item{edges}{which of a range's edges to include in the range}

\item{...}{ignored}
}
\value{
index vector or match matrix
}
\description{
test for \code{x} to lie between \code{low} and \code{high} (including the
margins), i.e., \eqn{low \leq x \leq high}
}
\details{
\code{between} is vectorized for \code{x} as well as for the limits. In the
case of vectorized limits \code{low} and \code{high} have to have the same
length and the ranges are defined by corresponding pairs from \code{low} and
\code{high} (those ranges are allowed to overlap).
}
\examples{
between(c(1,pi), 0:3,2:5+.5, TRUE)
between(c(1,pi), 0:3,2:5+.5, FALSE)
}
\author{
Benno Pütz \email{puetz@psych.mpg.de}
}
