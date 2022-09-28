#' Construct elliptical filters for matched CLIP-seq experiments
#'
#' This function constructs an elliptical filter centered on the origin for two
#' libraries based on the per nucleotide read counts.
#' The filter is defined as \code{((exp)^2)/a + ((ctrl)^2)/b = 1}, where
#' \code{a} and \code{b} are the means of each coverage vector plus \code{sdn}
#' standard deviations. Positions within the ellipse return \code{FALSE}, while
#' those outside the ellipse return \code{TRUE}. Will return an IRanges rle
#' vector if the input vectors are rle vectors.
#'
#' @param exp A (possibly IRanges rle) vector containing read-counts per
#' nucleotide in a crosslinked sample
#' @param ctrl A (possibly IRanges rle) vector containing read-counts per
#' nucleotide in a non-crosslinked control
#' @param sdn The number of standard deviations to use in constructing the
#' elliptical filter
#'
#' @return Returns a (possibly IRanges rle) elliptical filter mask.
#'
#'
#' @examples
#'
#' @export
#'
#'

filter_elliptical <- function(exp, ctrl, sdn=8){
  a <- (mean(exp[exp > 0]) + sdn*sd(exp[exp > 0]))^2
  b <- (mean(ctrl[ctrl > 0]) + sdn*sd(ctrl[ctrl > 0]))^2
  return(((exp)^2)/a + ((ctrl)^2)/b > 1)
}
