#' Find separation point between two clusters in a 1D vector
#'
#'
#'
#' @param data A 1-dimensional vector
#' @param plot Plot histogram?
#'
#' @return Returns value for separation point between clusters
#'
#' @details This function uses the \code{\link[Ckmeans.1d.dp]{Ckmeans.1d.dp}} package to
#' compute an optimal 2-component kmeans clustering on the given data, returning
#' the lowest value within the cluster with the highest mean to be used as a
#' separator in CLIP normalization.
#'
#'
#' @examples
#'
#' @seealso \code{\link{norclip}}, \code{\link{clipScaleFactors}},
#' \code{\link{Ckmeans.1d.dp}}
#'
#' @export
#' @import Ckmeans.1d.dp
#'
find_separator_kmeans <- function(data, plot=TRUE){
  clustering <- Ckmeans.1d.dp(data, 2)
  foreground <- which(clustering$centers == max(clustering$centers))
  fg_indices <- which(clustering$cluster == foreground)
  separator <- min(data[fg_indices])

  if(plot){
    hist(data, breaks="FD", main=paste("K-means clustering breakpoint"), col="lightgrey")
    abline(v=separator, col="red")
  }

  return(separator)
}
