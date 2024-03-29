% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{centering}
\alias{centering}
\alias{mc}
\alias{scale2}
\title{Utilities for centering/scaling}
\usage{
mc(x)

scale2(x, center = TRUE, scale = FALSE)
}
\arguments{
\item{x}{Numeric matrix}

\item{center, scale}{Passed to \code{scale(x, ...)}}
}
\description{
\code{mc(x)} is a short for \code{sweep(x, 2, colMeans(x))}

\code{scale2(x)} is essentially \code{scale(x)},
but returns 0 rather than \code{NaN} for constant columns
}
