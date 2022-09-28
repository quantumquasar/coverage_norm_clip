#' A function to generate diagnostic plots for CLIP-seq data.
#'
#' This function produces various diagnostic plots to help the user understand
#' if their data is suitable for \code{\link{norclip}} analysis, and if the
#' default filtering values are appropriate. These include correlations between
#' libraries, unfiltered and filtered 2D density plots comparing crosslinked
#' libraries to uncrosslinked control libraries, and 1D histograms of the
#' resulting log2 ratios. See vignette for advice on interpreting these plots.
#'
#'
#' @param wigs A list of IRanges rle coverage vectors.
#' @param data_table Vector metadata in data frame format, see package vignette
#' for details.
#' @param sdn Number of standard deviations for elliptical filtering, see
#' \code{\link{filter_elliptical}} for details.
#' @param colramp Color scheme for 2D density plots. Defaults to mimic the color
#' scheme of \code{\link{smoothScatter}}.
#'
#' @return No explicit return, will produce a correlation plot plus 3 plots per
#' experiment group: an unfiltered 2D density plot, an elliptically filtered 2D
#' density plot, and a density plot for the log2 ratio of experiment to control
#' library read counts.
#'
#' @examples
#'
#' @export
#'
#' @seealso \code{\link{norclip}}, \code{\link{loadData}},
#' \code{\link{filter_elliptical}}, \code{\link{clip_scat}}

runDiagnostics <- function(wigs, data_table, sdn=10,
                           colramp=rev(rainbow(10, end = 4/6))){
  message("Running CLIP diagnostics")
  colnames(data_table) <- c("identifier","type","direction","file")

  cor_mat <- plotCorrelations(wigs, data_table)

  uids <- as.vector(unique(data_table$identifier), mode="list")

  factors <- plyr::l_ply(uids, function(this_id){
    efi <- which(data_table$identifier == this_id & data_table$type == "E" &
                   data_table$direction =="F")
    eri <- which(data_table$identifier == this_id & data_table$type == "E" &
                   data_table$direction =="R")
    cfi <- which(data_table$identifier == this_id & data_table$type == "C" &
                   data_table$direction =="F")
    cri <- which(data_table$identifier == this_id & data_table$type == "C" &
                   data_table$direction =="R")

    erle <- c(wigs[[efi]], wigs[[eri]])
    crle <- c(wigs[[cfi]], wigs[[cri]])

    filt_0 <- which(erle > 0 | crle > 0)

    clip_scat(as.vector(erle[filt_0], mode="numeric"),
              as.vector(crle[filt_0], mode="numeric"), elliptical=0,
              main=paste("Unfiltered 2D frequency plot for", this_id))

    filt <- filter_elliptical(erle, crle, sdn=sdn)

    evec <- as.vector(erle[filt], mode="numeric")
    cvec <- as.vector(crle[filt], mode="numeric")
    rm(erle, crle)

    clip_scat(evec, cvec, elliptical = 0,
              main=paste("Filtered 2D frequency plot for", this_id))

    ratio <- log2(evec/cvec)

    hist(ratio, breaks="FD", freq=FALSE,
         main=paste("Ratio density for ", this_id, sep= " "), col="lightgrey")
  }, .progress="text")

}
