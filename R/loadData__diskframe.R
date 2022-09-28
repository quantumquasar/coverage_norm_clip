#' Load CLIP-seq data for analysis with \code{norclip}
#'
#' @param data_table A data frame describing the experiments.
#'
#' @return Returns a list of IRanges rle coverage vectors.
#'
#' @details This function loads data to be used for analysis in various
#' functions provided by the \code{norclip} package. The format of the input
#' data frame should be \code{identifier}, \code{type}, \code{direction},
#' \code{file}. See the vignette for details.
#'
#' @examples
#'
#' @export

loadData_df <- function(data_table, replicons){
  stopifnot(is.data.frame(data_table), ncol(data_table) == 4)
  colnames(data_table) <- c("identifier","type","direction","file")

  print("-> in loadData")
  #sanity check
  plyr::l_ply(unique(data_table$identifier), function(x){
    tmp <- data_table[data_table$identifier == x,]
    if(length(tmp[, "identifier"]) != 4){
      stop("Number of entries for identifer ",x," not equal to four.")
    }
    if(length(tmp[tmp$type == "E" & tmp$direction == "F",1]) != 1){
      stop("Incorrect number of forward experimental files for ", x)
    }
    if(length(tmp[tmp$type == "E" & tmp$direction == "R",1]) != 1){
      stop("Incorrect number of reverse experimental files for ", x)
    }
    if(length(tmp[tmp$type == "C" & tmp$direction == "F",1]) != 1){
      stop("Incorrect number of forward control files for ", x)
    }
    if(length(tmp[tmp$type == "C" & tmp$direction == "R",1]) != 1){
      stop("Incorrect number of reverse control files for ", x)
    }
  })

  message("Assuming ", length(unique(data_table$identifier)), " experiments")

  message("Reading wiggle files")

  start_loaddata <- Sys.time()


  #llply returns the filenames as a list
  # wigs <- plyr::llply(data_table$file, function(x){
  #   this_wig <- import(as.character(x))
  #
  #
  #
  #   this_wig <- abs(coverage(this_wig, weight="score"))
  #
  #
  #   return(this_wig)
  # }, .progress="text")



  #######
  ###using disk.frame

  #define df_path, temp dirs
  df_path <- getwd()


  wigs <- plyr::llply(data_table$file, function(x){
    wigs.df <- csv_to_disk.frame(x, fill=T, header = F,
                                 #keep.rownames=T,
                                 outdir = paste0(df_path,"tmp_all"))

    #print(format(object.size(wigs.df),"Kb"))

    wigs.df.2 <- wigs.df %>%
      select(V1, V2, V3)%>%
      collect()%>%
      mutate(id=row_number())%>%
      as.disk.frame(outdir = paste0(df_path, "tmp2"), overwrite = T)

    #print(paste0("df.2", format(object.size(wigs.df.2), "Kb")))

    start_locs <- as.numeric(wigs.df.2[V1=="variableStep", id])
    #start_locs

    chr_names <- wigs.df.2[V1=="variableStep", V2]
    #chr_names

    list_rle <- NULL
    for(y in 1:length(start_locs)){
      #is it possible to use llply, function instead of the for?

      ##calc RLEs for each chromosome
      if(!(y==length(start_locs))){

        pos <- wigs.df.2 %>%
          select(V1, V2, id) %>%
          collect()%>%
          slice((start_locs[y]+1):(start_locs[y+1]-1))

        rl3 <- IRanges(as.numeric(pos[, V1]),as.numeric(pos[, V1]))
        scores <- as.numeric(pos$V2)
        cov <- abs(coverage(rl3, weight = scores))
        if(isEmpty(list_rle)){
          list_rle <- RleList(cov)
        }
        else{
          list_rle <- append(list_rle, RleList(cov))
        }

      }
      else{
        pos <- wigs.df.2 %>%
          select(V1, V2, id) %>%
          collect()%>%
          slice(start_locs[y]+1:nrow(wigs.df.2))

        rl3 <- IRanges(as.numeric(pos[, V1]),as.numeric(pos[, V1]))
        scores <- as.numeric(pos$V2)
        cov <- abs(coverage(rl3, weight = scores))
        if(isEmpty(list_rle)){
          list_rle <- RleList(cov)
        }
        else{
          list_rle <- append(list_rle, RleList(cov))
        }
      }
      print(y)
    }

    names(list_rle)<-chr_names
    #list_rle

    return(list_rle)
  }, .progress="text")



  #####
  #check replicons identical, get lengths
  if(missing(replicons)){
    #cool, just comparing replicon_tags in the RLElist and
    #assigning names to replicons
    replicons <- names(wigs[[1]])
  }
  #print(replicons)
  print("->still in loadData1")


  #check if replicons extracted from expermient 1 are the same for all
  #the experiments
  plyr::l_ply(wigs, function(x){
    if(length(setdiff(replicons, names(x))) > 0){
      stop("Replicons differ between experiments")
    }
  })
  print("->still in loadData2")


  #implement sampling strategy for large genomes?

  ##this chunk assigns the maximum length of a replicon (RLElist element) to max_len
  #in this case, 4 replicons

  max_len <- as.list(rep(0, length(replicons)))
  max_len <- setNames(max_len, replicons)
  for (name in replicons){
    for (i in 1:length(wigs)){
      if(length(wigs[[i]][[name]]) > max_len[[name]]){
        max_len[[name]] <- length(wigs[[i]][[name]])
      }
    }
  }

  print("->still in loadData3")

  #adjust lengths
  for (name in replicons) {
    for (i in 1:length(wigs)){
      diff <- max_len[[name]] - length(wigs[[i]][[name]])
      if(diff > 0){
        wigs[[i]][[name]] <- c(wigs[[i]][[name]], rep(0, times=diff))
      }
    }
  }


  #convert all replicons to single RLE
  wigs <- plyr::llply(wigs, function(x){
    return(unlist(x[order(replicons)]))
  })
  print("->still in loadData4")

  return(wigs)
}

