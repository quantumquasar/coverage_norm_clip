% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotCorrelations.R
\name{plotCorrelations}
\alias{plotCorrelations}
\title{Plot clustered correlation heatmap for a list of RLE coverage plots}
\usage{
plotCorrelations(
  wigs,
  data_table,
  method = "pearson",
  log = F,
  heatscale = c(low = "darkblue", high = "lightblue")
)
}
\arguments{
\item{wigs}{A list of IRanges rle coverage vectors.}

\item{data_table}{Vector metadata in data frame format, see package vignette
for details.}

\item{method}{Method for cor function, defaults to "pearson"}

\item{log}{Log transform read counts?}

\item{heatscale}{Color range to use for heatscale.}
}
\value{
Returns a correlation matrix.
}
\description{
Plot clustered correlation heatmap for a list of RLE coverage plots
}
\details{
This function calculates correlation values between a set of
(possibly IRanges RLE) coverage plots, plotting them using the
ggplot2. Additionally, the correlation matrix is
returned.
}
\examples{

}
\seealso{
\code{\link{runDiagnostics}}
}
