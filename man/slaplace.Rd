% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/margins.R
\name{slaplace}
\alias{slaplace}
\alias{dslaplace}
\alias{pslaplace}
\alias{qslaplace}
\alias{rslaplace}
\title{Skew Laplace distribution}
\usage{
dslaplace(x, mu = 0.05, scale = 1, gamma = 1, log = FALSE)

pslaplace(q, mu = 0.05, scale = 1, gamma = 1)

qslaplace(p, mu = 0.05, scale = 1, gamma = 1)

rslaplace(n, mu = 0.05, scale = 1, gamma = 1)
}
\arguments{
\item{x}{vector of values}

\item{mu}{location parameter}

\item{scale}{scale parameter}

\item{gamma}{parameter}

\item{log}{flag for log density}

\item{q}{vector of quantiles}

\item{p}{vector of probabilities}

\item{n}{number of observations}
}
\description{
Skew Laplace distribution
}
