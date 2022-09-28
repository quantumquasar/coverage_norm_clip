#' Find local minima in the empirical distribution of a vector.
#'
#'
#' Heuristic function to find local minima in the empirical distribution of a
#' vector of values within a certain range. Usually used with the output of
#' \code{\link{find_maxima}} to identify a local minimum. Minima are identified
#' as appropriate change points in a density curve fitted to the input data, and
#' are reported in order of increasing density in the case of multiple minima.
#' Will optionally plot the density and predicted minima over a histogram of the
#' input data.
#'
#'
#' @param data A vector containing numeric data points.
#' @param n Number of local minima to find in the empirical distribution.
#' @param range Range in data within which to find minima.
#' @param plot Plot density and found minima?
#'
#' @return A vector containing n local minima ordered by
#' decreasing density.
#'
#' @examples
#'
#' x <- c(rnorm(100), rnorm(100, mean=4))
#' maxima <- find_maxima(x)
#' find_minima_in_range(x, range=c(min(maxima), max(maxima)), plot=TRUE)
#'
#' @seealso \code{\link{find_maxmima}}
#'
#' @export
#'
#'

find_minima_in_range <- function(data, n=1,
                                 range=c(min(data), max(data)), plot=FALSE, density_est="empirical", adjust=2, components=2){
  minima <- NULL

  density <- density(data, adjust=adjust)
  if(density_est == "model"){
    mod = mclust::densityMclust(data, G=2)
    summary(mod)
    density <- as.data.frame(cbind(x=density$x, y=predict(mod, density$x)))
  }

  density_range <- density$x > range[1] & density$x < range[2]

  y <- density$y[density_range]
  spacer <- length(which(density$x <= range[1]))

  for ( i in 2:(length(y)-1) ){
    if ( (y[i] < y[i-1]) & (y[i] < y[i+1]) ) {
      minima <- c(minima,i + spacer)
    }
  }
  if ( length(minima) == 0 ) {
    stop('Distribution appears monotonic')
  }

  if ( length(minima) < n ) {
    stop(paste('Only',length(minima), "minima found!", sep=" "))
  }

  positions <- density$x[minima]
  densities <- density$y[minima]

  order <- order(densities, decreasing=FALSE)

  selected_positions <- positions[order[1:n]]

  if(plot){
  #  if(density_est=="empirical"){
      hist(data, breaks="FD", freq=FALSE,
           main=paste("First",n, "local minima", sep= " "), col="lightgrey")
      lines(density, lwd=2)
      plyr::l_ply(selected_positions, function(x){abline(v=x, col="red", lwd=2)})
  #  }
  #  else if(density_est=="model"){
   #   plot(mod, what="density", data=data, breaks="FD",
  #         main=paste("First",n, "local minima", sep= " "), col="lightgrey")
  #    plyr::l_ply(selected_positions, function(x){abline(v=x, col="red", lwd=2)})
   # }
  }

  return(selected_positions)
}
