% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadData__fread.R
\name{loadData_fread}
\alias{loadData_fread}
\title{Load CLIP-seq data for analysis with \code{norclip}}
\usage{
loadData_fread(data_table, replicons)
}
\arguments{
\item{data_table}{A data frame describing the experiments.}
}
\value{
Returns a list of IRanges rle coverage vectors.
}
\description{
Load CLIP-seq data for analysis with \code{norclip}
}
\details{
This function loads data to be used for analysis in various
functions provided by the \code{norclip} package. The format of the input
data frame should be \code{identifier}, \code{type}, \code{direction},
\code{file}. See the vignette for details.
}
\examples{

}
