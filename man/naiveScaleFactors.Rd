% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/naiveScaleFactors.R
\name{naiveScaleFactors}
\alias{naiveScaleFactors}
\title{Calculate naive scale factors for RNA-seq data}
\usage{
naiveScaleFactors(wigs, data_table, sdn = 5, bg_cut = 5, plot = F)
}
\arguments{
\item{wigs}{A list of IRanges rle coverage vectors.}

\item{data_table}{Vector metadata in data frame format, see package vignette
for details.}

\item{sdn}{Number of standard deviations for elliptical filtering, see
\code{\link{filter_elliptical}} for details.}

\item{bg_cut}{Cut-off for library normalization.}

\item{plot}{Produce diagnostic plots?}
}
\value{
Returns an array of scale factors for the provided libraries.
}
\description{
Calculate naive scale factors for RNA-seq data
}
\details{
This function computes scale factors using a 'naive' method, i.e.
not considering the different sizes of experiment and control RNA pools.
This method was adapted from RNA-seq analysis, propsed by Anders and Huber,
Genome Biology, 2010.
}
\examples{

}
\seealso{
\code{\link{gm_scale_factors}}
}
