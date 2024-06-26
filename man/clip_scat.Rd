% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clip_scat.R
\name{clip_scat}
\alias{clip_scat}
\title{A function for producing density scatterplots for CLIP-seq data}
\usage{
clip_scat(
  exp,
  ctrl,
  sf_exp = 1,
  sf_ctrl = 1,
  xyline = F,
  exp_lab = "crosslinked",
  ctrl_lab = "non-crosslinked",
  main = "2D density plot",
  elliptical = 2,
  nbin = 100,
  colramp = rev(rainbow(10, end = 4/6))
)
}
\arguments{
\item{exp}{A vector containing read-counts per nucleotide in a crosslinked
sample}

\item{ctrl}{A vector containing read-counts per nucleotide in a
non-crosslinked control}

\item{sf_exp}{A scale factor to scale the crosslinked sample by.}

\item{sf_ctrl}{A scale factor to scale the non-crosslinked sample by.}

\item{xyline}{Draw the line \code{x = y} in red?}

\item{exp_lab}{A label describing the experimental condition.}

\item{ctrl_lab}{A label describing the control condition.}

\item{main}{A plot title.}

\item{elliptical}{The number of standard deviations to use in constructing an
elliptical filter.}

\item{nbin}{Number of bins to use for calculating hex density.}

\item{colramp}{Color scheme; defaults to imitating the
\code{\link{smoothScatter}} colors.}
}
\description{
Draws a 2D density plot for per-nucleotide read counts using the
\code{\link{ggplot}} function \code{\link{geom_hex}}.
}
\examples{

}
\seealso{
\code{\link{filter_elliptical}} \code{\link{geom_hex}}
}
