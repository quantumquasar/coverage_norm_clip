#' Find local maxima in the empirical distribution of a vector of values
#'
#' Heuristic function to find local maxima in the empirical distribution of a
#' vector of values within a certain range. Usually used to identify a range
#' for \code{\link{find_minima_in_range}} in order to identify a local minimum.
#' Minima are identified as appropriate change points in a density curve fitted
#' to the input data, and are reported in order of decreasing density in the
#' case of multiple maxima.  Will optionally plot the density and predicted
#' maxima over a histogram of the input data.
#'
#' @param data A vector containing numeric data points.
#' @param n Number of local maxima to find in the empirical distribution.
#' @param plot Plot density and maxima?
#'
#' @return A vector containing n local maxima ordered by
#' decreasing density.
#'
#' @examples
#'
#' x <- c(rnorm(100), rnorm(100, mean=4))
#' find_maxima(x, plot=TRUE)
#'
#' @seealso \code{\link{find_minima_in_range}}
#'
#' @export
#'

find_maxima <- function(data, n=2, plot=FALSE, density_est="model", adjust=2, components=2){
  maxima <- NULL

  density <- density(data, adjust=adjust)
  if(density_est == "model"){
    mod = mclust::densityMclust(data, G=components)
    density <- as.data.frame(cbind(x=density$x, y=predict(mod, density$x)))
  }
  #density <- density(data)
  y <- density$y
  for ( i in 2:(length(y)-1) ){
    if ( (y[i] > y[i-1]) & (y[i] > y[i+1]) ) {
      maxima <- c(maxima,i)
    }
  }
  if ( length(maxima) == 0 ) {
    stop('Distribution appears monotonic')
  }

  if ( length(maxima) < n ) {
    stop(paste('Only',length(maxima), "maxima found!", sep=" "))
  }

  positions <- density$x[maxima]
  densities <- density$y[maxima]

  order <- order(densities, decreasing=TRUE)

  selected_positions <- positions[order[1:n]]

  if(plot){
      hist(data, breaks="FD", freq=FALSE,
         main=paste("First",n, "local maxima", sep= " "), col="lightgrey")
      lines(density, lwd=2)
      plyr::l_ply(selected_positions, function(x){abline(v=x, col="red", lwd=2)})
    #else if(density_est=="model"){
     # plot(mod, what="density", data=data, breaks="FD",
      #     main=paste("First",n, "local maxima", sep= " "), col="lightgrey")
      #plyr::l_ply(selected_positions, function(x){abline(v=x, col="red", lwd=2)})
    #}
  }

  return(selected_positions)
}
