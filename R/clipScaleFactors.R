#' Calculate scale factors for CLIP-seq data
#'
#'
#'
#' @param wigs A list of IRanges rle coverage vectors.
#' @param data_table Vector metadata in data frame format, see package vignette
#' for details.
#' @param sdn Number of standard deviations for elliptical filtering, see
#' \code{\link{filter_elliptical}} for details.
#' @param crossnormalize Crossnormalize the control libraries?
#' @param plot Produce diagnostic plots?
#' @param bg_cut Cut-off for control library normalization.
#'
#' @return Returns an array of scale factors for the provided libraries.
#'
#' @details This function computes scale factors for CLIP-seq data, under the
#' assumption that it follows a bi- or multi-modal distribution. For each
#' matched experiment, low count reads are first filtered using an elliptical
#' filter based on the read count standard deviation. Log2 ratios of read counts
#' per nucleotide are calculated, and a local minima is identified between the
#' two highest density maxima. Positions with a log2 ratio less than this minima
#' are then used to fit a scale factor for each experiment. Optionally using the
#' \code{crossnormalize} parameter, all libraries can be scaled based on scale
#' factors computed between control libraries to produce consistent values for
#' hypothesis testing. Note positions with less than \code{bg_cut} reads in
#' any library will be excluded.
#'
#'
#' @examples
#'
#' @seealso \code{\link{norclip}}, \code{\link{loadData}},
#' \code{\link{gm_scale_factors}}, \code{\link{filter_elliptical}}
#' \code{\link{clip_scat}}
#'
#' @import IRanges
#' @export
#'

clipScaleFactors <- function(wigs, data_table, sdn=5, crossnormalize=T,
                             plot=T, bg_cut=5, breakpoint="kmeans", adjust=2){
  colnames(data_table) <- c("identifier","type","direction","file")
  uids <- as.vector(unique(data_table$identifier), mode="list")

  message("Calculating size factors")

  fg_sfs <- plyr::laply(uids, function(this_id){
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

    filt <- filter_elliptical(erle, crle, sdn=sdn)

    evec <- as.vector(erle[filt], mode="numeric")
    cvec <- as.vector(crle[filt], mode="numeric")
    rm(erle, crle)

    ratio <- log2((evec+1) / (cvec+1))

    #if(density_est == "empirical"){
    #  density <- density(ratio)
    #}
    #else if(density_est == "model"){
    #  mod = densityMclust(ratio, G=2)
    #  density <- as.data.frame(cbind(x=as.numeric(ratio), y=as.numeric(mod$density)))
    #}
    if((breakpoint == "empirical") | breakpoint == "model"){
      maxima <- find_maxima(ratio, plot=plot, density_est=breakpoint, adjust=adjust)

      separator <- find_minima_in_range(ratio, range=sort(maxima), plot=plot, density_est=breakpoint, adjust=adjust)
    } else if(breakpoint == "kmeans"){
      separator <- find_separator_kmeans(ratio, plot=plot)
    }
    scale_indices <- which(ratio < separator)

    message(paste(this_id, ":", length(scale_indices),
                  " background positions used for normalization", sep=""))

    sf <- median(evec[scale_indices] / cvec[scale_indices])

    if(plot){
      nf <- median(evec / cvec)
      #plot(sort(evec[scale_indices] / cvec[scale_indices]),
      #     ylab="XL+ to XL- ratio", xlab="Sort Index", pch=20)
      #abline(h=nf, col="red")
      #abline(h=sf, col="blue")
      #mtext(paste("naive scale factor:", nf), side=3, adj=0, line=0, col="red")
      #mtext(paste("norclip scale factor:",sf), side=3, adj=0, line=1,
      #      col="blue")
      #mtext(this_id, side=3, adj=0, line=2)

      plot(sort(evec / cvec),
           ylab="XL+ to XL- ratio", xlab="Sort Index", ylim=c(0,5*nf), pch=20)
      abline(h=nf, col="red")
      abline(h=sf, col="blue")
      mtext(paste("naive scale factor:", nf), side=3, adj=0, line=0, col="red")
      mtext(paste("norclip scale factor:",sf), side=3, adj=0, line=1,
            col="blue")
      mtext(paste(this_id, "- no filter"), side=3, adj=0, line=2)
    }

    return(sf)
  })

  names(fg_sfs) <- unlist(uids)
  #return(fg_sfs)



  if(crossnormalize){
    bg_vecs <- plyr::llply(uids, function(this_id){
      cfi <- which(data_table$identifier == this_id & data_table$type == "C" &
                     data_table$direction =="F")
      cri <- which(data_table$identifier == this_id & data_table$type == "C" &
                     data_table$direction =="R")
      crle <- c(wigs[[cfi]], wigs[[cri]])
      return(crle)
    })

    #return(bg_vecs)

    nz_pos <- plyr::llply(bg_vecs, function(this_rle){
      return(which(this_rle > bg_cut))
    })

    nz <- Reduce(intersect, nz_pos)
    rm(nz_pos)
    message(paste(length(nz),
                  " background positions used for crossnormalization", sep=""))

    nz_ar <- plyr::laply(bg_vecs, function(this_rle){
      return(as.vector(this_rle[nz], mode="numeric"))
    })

    bg_sfs <- gm_scale_factors(t(nz_ar))
    names(bg_sfs) <- uids
  } else {
    bg_sfs <- rep(1, times=length(uids))
  }

  names(bg_sfs) <- unlist(uids)

  sfs <- plyr::aaply(seq(from=1,to=length(unlist(uids))), 1, function(this_ind){
    return(cbind(bg_sfs[this_ind] * fg_sfs[this_ind], bg_sfs[this_ind]))
  })
  rownames(sfs) <- unlist(uids)
  colnames(sfs) <- c("E", "C")

  if(plot){
    plyr::l_ply(uids, function(this_id){
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

      message(paste(this_id,sfs[this_id,"E"], sfs[this_id,"C"]))
      clip_scat(erle, crle, sf_exp=sfs[this_id,"E"], sf_ctrl=sfs[this_id, "C"],
                xyline=T, elliptical=sdn,
                main=paste("2D frequency plot for", this_id))
    })
  }

  return(sfs)
}
