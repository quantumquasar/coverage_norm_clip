#' Plot clustered correlation heatmap for a list of RLE coverage plots
#'
#'
#'
#' @param wigs A list of IRanges rle coverage vectors.
#' @param data_table Vector metadata in data frame format, see package vignette
#' for details.
#' @param method Method for cor function, defaults to "pearson"
#' @param log Log transform read counts?
#' @param heatscale Color range to use for heatscale.
#'
#' @return Returns a correlation matrix.
#'
#' @details This function calculates correlation values between a set of
#' (possibly IRanges RLE) coverage plots, plotting them using the
#' ggplot2. Additionally, the correlation matrix is
#' returned.
#'
#'
#' @examples
#'
#' @seealso \code{\link{runDiagnostics}}
#'
#' @export


plotCorrelations <- function(wigs, data_table, method="pearson", log=F,
                             heatscale= c(low='darkblue',high='lightblue')){
  colnames(data_table) <- c("identifier","type","direction","file")
  nz_pos <- plyr::llply(wigs, function(this_rle){return(which(this_rle != 0))})
  nz <- Reduce(union, nz_pos)

  vec_list <- plyr::llply(wigs, function(this_rle){as.vector(this_rle[nz],
                                                             mode="numeric")})
  nz_mat <- do.call(rbind, vec_list)
  these_names <- c()
  for(index in 1:length(data_table[,1])){
    these_names <- c(these_names,
                     paste(data_table$identifier[index], data_table$type[index],
                           data_table$direction[index], sep="."))
  }
  if(log){
    cormat <- cor(t(log(nz_mat, base=2)), method="pearson")
  } else {
    cormat <- cor(t(nz_mat), method="pearson")
  }
  rownames(cormat) <- these_names
  colnames(cormat) <- these_names


  # Adapted from a snippet found at:
  # https://www.r-bloggers.com/ggheat-a-ggplot2-style-heatmap-function/

  m <- cormat[hclust(dist(cormat))$order ,hclust(dist(t(cormat)))$order]
  rows <- dim(m)[1]
  cols <- dim(m)[2]
  melt.m <- cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows),
               reshape2::melt(m))
  g <- ggplot2::ggplot(data=melt.m)

  g <- g + ggplot2::geom_rect(ggplot2::aes(xmin=colInd-1,xmax=colInd,
                                       ymin=rowInd-1,ymax=rowInd, fill=value)
                          ,colour='white')

  g <- g + ggplot2::scale_x_continuous(breaks=(1:cols)-0.5, labels=colnames(m))
  g <- g + ggplot2::scale_y_continuous(breaks=(1:rows)-0.5, labels=rownames(m))

  g <- g + ggplot2::theme(panel.grid.minor=ggplot2::element_line(colour=NA),
                       panel.grid.major=ggplot2::element_line(colour=NA),
                       panel.background=ggplot2::element_rect(fill=NA,
                                                              colour=NA))

  print(g)
  return(cormat)
}
