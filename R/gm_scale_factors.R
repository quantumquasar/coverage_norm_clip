#' A function for determining RNA-seq scale factors using the geometric mean
#' across samples.
#'
#' This function implements an RNA-seq normalization procedure based on the
#' median read ratio between a sequencing library and a pseudoreference
#' constructed as the geometric mean across features, as introduced by Anders
#' and Huber, Genome Biology, 2010.
#'
#'
#' @param matrix A matrix containing columns corresponding to samples, and rows
#' to per feature read counts
#' @param locfunc Location function for normalization; defaults to median.
#'
#' @return An array of scale factors.
#'
#' @examples
#'
#' #trivial example
#' counts <- sample(2:5000, size=200, replace=T)
#' count_mat <- matrix(c(counts, 2 * counts, 3 * counts), ncol=3)
#' sfs <- gm_scale_factors(count_mat)
#' sum(count_mat[,1] / sfs[1] - count_mat[,2] / sfs[2])
#' # should equal 0
#'
#' d_1 <- sample(-2:2, size=200, replace=T)
#' d_2 <- sample(-2:2, size=200, replace=T)
#' count_mat <- matrix(c(counts, 2 * (counts + d_1), 3 *  (counts+d_2)), ncol=3)
#' sfs_d <- gm_scale_factors(count_mat)
#' sfs - sfs_d
#' # These values should be small, showing scale factors should robust to small
#' # deviations.
#'
#'
#' @export
#'
#'

gm_scale_factors <- function(matrix, locfunc=median){
  log_mean <- rowMeans(log(matrix))
  nz <- which(log_mean != 0 & !is.infinite(log_mean))
  SF <- plyr::aaply(matrix[nz,],2,function(x){exp(locfunc(log(x) - log_mean[nz]))})

  return(SF)
}
