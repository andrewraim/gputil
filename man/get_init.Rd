% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init.R
\name{get_init}
\alias{get_init}
\title{Gibbs initial values object}
\usage{
get_init(n, muX = NULL, sigma2 = NULL)
}
\arguments{
\item{n}{Number of observations.}

\item{muX}{A vector of length \code{n}.}

\item{sigma2}{A scalar}
}
\value{
An object of class \code{my_init}.
}
\description{
Gibbs initial values object
}
