#' A function for calculating the geometric mean of a vector
#'
#'
#' @param x A vector containing
#' @param na.rm Remove undefined values?
#' @param zero.propagate Propagate zeros?
#'
#' @return gm_mean The geometric mean of a vector.
#'
#' @examples
#'
#' @export
#'
#'


gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}
