% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistical_tools.R
\name{pracf}
\alias{pracf}
\title{Partial Rank Autocorrelation Function}
\usage{
pracf(
  data,
  method = "kendall",
  nlags = floor(5 * log10(length(data))),
  plot = TRUE,
  main = paste("Series", deparse(substitute(data)), sep = " "),
  xlab = "Lag",
  ylab = "PACF",
  lwd = 2
)
}
\arguments{
\item{data}{a vector of data to which the autocorrelation function is to be computed}

\item{method}{one of three ways to compute the acf: "kendall", "spearman" or "pearson"}

\item{nlags}{the number of lags considered}

\item{plot}{logical specifing whether to plot the acf or not}

\item{main}{title of the plot}

\item{xlab}{x-axis label of the plot}

\item{ylab}{y-axis label of the plot}

\item{lwd}{width of lines}

\item{exact}{logical indicating whether the rpacf should be exact or fast}
}
\value{
a list containing the partial autocorrelation values and pvalues
}
\description{
Partial Rank Autocorrelation Function
}
\examples{
pracf(runif(1000))
y <- sim(dvinecopula("gauss", list(0.5)))
pracf(y)
}
