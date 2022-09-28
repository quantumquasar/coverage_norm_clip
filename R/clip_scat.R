#' A function for producing density scatterplots for CLIP-seq data
#'
#' Draws a 2D density plot for per-nucleotide read counts using the
#' \code{\link{ggplot}} function \code{\link{geom_hex}}.
#'
#' @param exp A vector containing read-counts per nucleotide in a crosslinked
#' sample
#' @param ctrl A vector containing read-counts per nucleotide in a
#' non-crosslinked control
#' @param sf_exp A scale factor to scale the crosslinked sample by.
#' @param sf_ctrl A scale factor to scale the non-crosslinked sample by.
#' @param xyline Draw the line \code{x = y} in red?
#' @param exp_lab A label describing the experimental condition.
#' @param ctrl_lab A label describing the control condition.
#' @param main A plot title.
#' @param elliptical The number of standard deviations to use in constructing an
#' elliptical filter.
#' @param nbin Number of bins to use for calculating hex density.
#' @param colramp Color scheme; defaults to imitating the
#' \code{\link{smoothScatter}} colors.
#'
#' @examples
#'
#' @export
#'
#' @seealso \code{\link{filter_elliptical}} \code{\link{geom_hex}}
#'

clip_scat <- function(exp, ctrl, sf_exp=1, sf_ctrl=1, xyline=F,
                      exp_lab="crosslinked", ctrl_lab="non-crosslinked",
                      main="2D density plot", elliptical=2, nbin=100,
                      colramp=rev(rainbow(10, end = 4/6))){
  if(elliptical > 0){
    filt <- filter_elliptical(exp, ctrl, elliptical)
    masked_exp <- as.vector(exp[filt], mode="numeric")
    masked_ctrl <- as.vector(ctrl[filt], mode="numeric")
  }else{
    masked_exp <- as.vector(exp, mode="numeric")
    masked_ctrl <- as.vector(ctrl, mode="numeric")
  }
  masked_exp <- masked_exp / sf_exp
  masked_ctrl <- masked_ctrl / sf_ctrl

  this_dat <- data.frame(exp = masked_exp, ctrl = masked_ctrl)

  p <- ggplot2::ggplot(this_dat, ggplot2::aes(x=exp, y=ctrl))
  p <- p + ggplot2::geom_hex(ggplot2::aes(x=exp,y=ctrl), bins=100)
  p <- p + ggplot2::xlab("+XL") + ggplot2::ylab ("-XL")
  p <- p + ggplot2::scale_fill_gradientn("", colours = colramp)
  p <- p + ggplot2::theme_bw()
  p <- p + ggplot2::ggtitle(main)
  if(xyline){
    p <- p + ggplot2::geom_abline(intercept=0, slope=1, colour="red")
  }
  print(p)
}
