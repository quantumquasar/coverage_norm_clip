#' Calculate naive scale factors for RNA-seq data
#'
#'
#' @param wigs A list of IRanges rle coverage vectors.
#' @param data_table Vector metadata in data frame format, see package vignette
#' for details.
#' @param sdn Number of standard deviations for elliptical filtering, see
#' \code{\link{filter_elliptical}} for details.
#' @param plot Produce diagnostic plots?
#' @param bg_cut Cut-off for library normalization.
#'
#' @return Returns an array of scale factors for the provided libraries.
#'
#' @details This function computes scale factors using a 'naive' method, i.e.
#' not considering the different sizes of experiment and control RNA pools.
#' This method was adapted from RNA-seq analysis, propsed by Anders and Huber,
#' Genome Biology, 2010.
#'
#' @examples
#'
#' @seealso \code{\link{gm_scale_factors}}
#'
#' @export

naiveScaleFactors <- function(wigs, data_table, sdn=5, bg_cut=5, plot=F){
  colnames(data_table) <- c("identifier","type","direction","file")
  uids <- as.vector(unique(data_table$identifier), mode="list")

  vecs <- plyr::llply(uids, function(this_id){
    efi <- which(data_table$identifier == this_id & data_table$type == "E" &
                   data_table$direction =="F")
    eri <- which(data_table$identifier == this_id & data_table$type == "E" &
                   data_table$direction =="R")
    erle <- c(wigs[[efi]], wigs[[eri]])

    cfi <- which(data_table$identifier == this_id & data_table$type == "C" &
                   data_table$direction =="F")
    cri <- which(data_table$identifier == this_id & data_table$type == "C" &
                   data_table$direction =="R")
    crle <- c(wigs[[cfi]], wigs[[cri]])
    return(list(exp=erle,ctrl=crle))
  })

  names(vecs) <- unlist(uids)

  #return(vecs)

  nz_pos_exp <- plyr::llply(vecs, function(this_uid){
    return(which(this_uid$exp > bg_cut))
  })

  nz_exp <- Reduce(intersect, nz_pos_exp)
  rm(nz_pos_exp)

  nz_pos_ctrl <- plyr::llply(vecs, function(this_uid){
    return(which(this_uid$ctrl > bg_cut))
  })

  nz_ctrl <- Reduce(intersect, nz_pos_ctrl)
  rm(nz_pos_ctrl)

  nz <- union(nz_exp, nz_ctrl)

  message(paste(length(nz)," background positions used for naive normalization",
                sep=""))

  vecs <- unlist(vecs)

  nz_ar <- plyr::laply(vecs, function(this_rle){
    return(as.vector(this_rle[nz], mode="integer"))
  })

  sfs <- gm_scale_factors(t(nz_ar))
  sfs <- matrix(sfs, ncol=2, byrow=T)
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
                main=paste("Naive 2D density plot for", this_id))
    })
  }

  return(sfs)
}
