#' Perform norclip normalization of CLIP-seq libraries
#'
#'
#' @param data_table A description table
#' @param crossnormalize Perform cross-normalization of replicates?
#' @param sdn Number of standard deviations for elliptical filtering
#' @param diagnostics Run diagnostics prior to normalization?
#' @param bg_cut Number of reads required in background library for a position
#' to be used in normalization
#' @param naive Additionally run "naive" normalization?
#'
#' @examples
#'
#' @export
#'
#'
#'

#data_table <-read.table(file, stringsAsFactors = F)
# colnames(data_table) <- c("identifier","type","direction","file")

norclip <- function(data_table, crossnormalize=T, sdn=5, diagnostics=T,
                    bg_cut=5, naive=T, export_wigs=F, adjust=2, breakpoint="kmeans", replicons){

  start_t <- Sys.time()
  #wigs <- loadData(data_table, replicons)

  #using fread
  wigs <- loadData_fread(data_table, replicons)
  end_t <- Sys.time()
  elapsed <- end_t - start_t
  cat('time to load: ', elapsed)


  if(diagnostics){
    runDiagnostics(wigs, data_table, sdn=sdn)
  }

  factors <- clipScaleFactors(wigs, data_table, sdn, crossnormalize,
                              plot=diagnostics, bg_cut=bg_cut, breakpoint=breakpoint, adjust=adjust)

  if(naive){
    naive <- naiveScaleFactors(wigs, data_table, sdn=sdn, bg_cut=bg_cut,
                               plot=diagnostics)
  }

  return(factors)
}
