#' Import wig files into rle coverage vectors.
#'
#' Imports wiggle files for forward and reverse strands of a CLIP-seq
#' experiment. These are then converted to rle coverage vectors, and the
#' absolute value of the reverse strand is taken to guard against negative
#' values before concatenation.
#'
#'
#' @param forward_path Path to forward strand wig file.
#' @param reverse_path Path to reverse strand wig file.
#'
#' @return Return a concatenated IRanges rle coverage vector.
#'
#' @examples
#'
#' @seealso \code{\link{import}}, \code{\link{loadData}}
#'

wig2rle <- function(forward_path, reverse_path){
  x <- import(forward_path)
  x <- coverage(x, weight="score")

  y <- import(reverse_path)
  y <- coverage(y, weight="score")

  r <- unlist(c(x,  abs(y)), use.names=F)
  return(unlist(r))
}
