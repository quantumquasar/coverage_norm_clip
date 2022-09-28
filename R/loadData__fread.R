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

loadData_fread <- function(data_table, replicons){
  stopifnot(is.data.frame(data_table), ncol(data_table) == 4)
  colnames(data_table) <- c("identifier","type","direction","file")

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


  #####new bit of code####
  #using chunks

  ##if two or replicons are in one chunk:
  ###get line number, repeat
  add_chunks <- function(x, i, t){
    y <- tryCatch(
      expr    = {
        if (t==1){#for all chunks
          fread(x,nrows = 10000,skip = 10000*i,header = F, fill = T, data.table = T)}
        else{#for the final chunk, which may have differing num of rows
          fread(x, nrows = (chunks*10000-n_rows),skip = 10000*i,
                header = F, fill = T, data.table = T)
        }
      },
      warning = function(warn){
        #cat('Warning: ', warn$message, '\n\n');
        n_line <- as.numeric(gsub('Stopped early on line (\\d+)\\..*','\\1',warn$message))
        if (!is.na(n_line)) {
          cat('replicon on ', n_line,'\n')
          if (t==1){#for all chunks
            #until the replicon (-1)
            chunk_1 <- fread(x, skip = 10000*i, nrows = (n_line-10000*i-1),
                               header = F, fill = T, data.table = T)
            #including the replicon until the end/next replicon (which will be caught)
            chunk_2 <- fread(x, skip = (n_line-1), nrows = (10000*(i+1)-n_line+1),
                               header = F, fill = T, data.table = T)
            return(rbind(chunk_1, chunk_2, fill=T))
          }
          else{#for the final chunk
            chunk_1 <- fread(x, skip = 10000*i, nrows = (n_line-10000*i-1),
                               header = F, fill = T, data.table = T)
            chunk_2 <- fread(x, skip = (n_line-1), nrows = (chunks*10000-n_line+1),
                               header = F, fill = T, data.table = T)
            return(rbind(chunk_1, chunk_2, fill=T))
            ##problem being there can be two or more variablesteps in one chunk
            #which will be caught
          }
        }
      }
    )
    return(y)
  }

  #future: optimise
  #break down in functions
  calc_rle <- function(wig_chunk, ch){

    df <- as.data.frame(wig_chunk)
    df$id <- rownames(df)
    start <- as.numeric(df[df$V1=="variableStep" | df$V1=="fixedStep",]$id)

    #initialize coord
    if (coord==0 & (ch==0)){
      coord <<- as.numeric(df[df$id==start[1],]$V1)
    }

    #if there's only one replicon in the chunk
    if (length(start) == 1){
      chr_names <<- append(chr_names, df[start,]$V2)
      if (start == 1 | start == 2){#where is the replicon name?
        #if at the beginning, then start from there

        if (is.null(list_cov)){
          #list_cov stores the coverage for a replicon and is used
          #to concatenate further rles for the same replicon

          rl4 <- IRanges(as.numeric(df[(start+1):nrow(df), ]$V1),
                         as.numeric(df[(start+1):nrow(df), ]$V1))

          scores <- as.numeric(df[(start+1):nrow(df), ]$V2)

          #calculate coverage normally if the chunk contains a replicon
          cov <- abs(coverage(rl4, weight = scores))

          #concatenate coverage of this chunk to the previous
          #this way runLengths are merged, if runValues are common between
          #the end of the previous chunk and the beginning of the current chunk
          list_cov <<- cov
        }
        else{
          if(start == 2){
            #check if positions between chunks have a discontinuity
            #if so, create an RLE for the missing positions
            #assuming that all data in wig format start with a "track type..." line
            #   #this predicament won't arise with the first replicon
            #   #for subsequent replicons, if the replicon name is in line 2 of the chunk
            #   #we'll just check for line 1 and pos (coord) of the last line of the prev chunk

            if (!(ch==0) & !(df[df$id==1,]$V1 == (coord+1))){
              list_cov <<- c(list_cov,
                             abs(coverage(IRanges(c((coord+1):(as.numeric(df[df$id==1,]$V1)-1)),c((coord+1):(as.numeric(df[df$id==1,]$V1)-1))),
                                          shift = -(coord), weight = rep(0, (as.numeric(df[df$id==1,]$V1) - coord)-1))))
            }
            if (!(ch==0)){
              rl4 <- IRanges(as.numeric(df[1:(start-1), ]$V1),
                             as.numeric(df[1:(start-1), ]$V1))

              scores <- as.numeric(df[1:(start+1), ]$V2)

              cov <- abs(coverage(rl4, shift = -(as.numeric(df$V1[1])-1), weight = scores))
              list_cov <<- c(list_cov, cov)
            }
          }


          #create Rlelist element
          if(isEmpty(list_rle)){
            list_rle <<- RleList(list_cov)
          }
          else{
            list_rle <<- append(list_rle, RleList(list_cov))
          }
          list_cov <<- NULL
          rl4 <- IRanges(as.numeric(df[(start+1):nrow(df), ]$V1),
                         as.numeric(df[(start+1):nrow(df), ]$V1))

          scores <- as.numeric(df[(start+1):nrow(df), ]$V2)
          cov <- abs(coverage(rl4, weight = scores))

          list_cov <<- cov
        }
      }

      else{#if start is further down in the chunk
        ##future:
        ##may have to check if there's another line above
        ##variableStep chrom=...
        rl4 <- IRanges(as.numeric(df[1:(start-1), ]$V1),
                       as.numeric(df[1:(start-1), ]$V1))

        scores <- as.numeric(df[1:(start-1), ]$V2)
        cov <- abs(coverage(rl4, shift = -(as.numeric(df$V1[1])-1), weight = scores))


        if (!(ch==0) & !(df[df$id==1,]$V1 == (coord+1))){
          cov <- c(abs(coverage(IRanges(c((coord+1):(as.numeric(df[df$id==1,]$V1)-1)),c((coord+1):(as.numeric(df[df$id==1,]$V1)-1))),
                                shift = -(coord), weight = rep(0, (as.numeric(df[df$id==1,]$V1) - coord)-1))), cov)
        }

        list_cov <<- c(list_cov, cov)

        #update RLElist
        if(isEmpty(list_rle)){
          list_rle <<- RleList(list_cov)
        }
        else{
          list_rle <<- append(list_rle, RleList(list_cov))
        }
        list_cov <<- NULL

        #create RLE for new replicon
        #check if the replicon isnt' the last line of the chunk
        #chr_names <- append(chr_names, df[start[s],]$V2)

        if ((start) < nrow(df)){
          rl4 <- IRanges(as.numeric(df[(start+1):nrow(df), ]$V1),
                         as.numeric(df[(start+1):nrow(df), ]$V1))

          scores <- as.numeric(df[(start+1):nrow(df), ]$V2)
          cov <- abs(coverage(rl4, weight = scores))

          list_cov <<- cov
        }
      }
    }
    #if there are no replicons
    else if (length(start) == 0){
      rl4 <- IRanges(as.numeric(df[(1):nrow(df), ]$V1),
                     as.numeric(df[(1):nrow(df), ]$V1))

      scores <- as.numeric(df[(1):nrow(df), ]$V2)

      #the coverage computation needs to be adjusted for the starting index
      #coverage assumes all indices less than the starting index to be 0s since these
      #are not present
      #cov2 <- abs(coverage(rl4, shift = -(df$V1[1]-1), weight = scores))
      cov2 <- abs(coverage(rl4, shift = -(df$V1[1]-1), weight = scores))
      #if_else(is.null(list_cov), list_cov = cov2, list_cov = c(list_cov, cov2))

      if (!(ch==0) & !(df[df$id==1,]$V1 == (coord+1))){
        cov2 <- c(abs(coverage(IRanges(c((coord+1):(as.numeric(df[df$id==1,]$V1)-1)),c((coord+1):(as.numeric(df[df$id==1,]$V1)-1))),
                               shift = -(coord), weight = rep(0, (as.numeric(df[df$id==1,]$V1) - coord)-1))), cov2)
      }

      list_cov <<- c(list_cov, cov2)
    }

    #if there are more than one replicon
    else{
      for (s in 1:length(start)){
        repl <- df[df$id==start[s],]$V2
        #run a loop to repl[id]+1 to repl[id==start[s+1]]-1
        #to cover the indices and values of this replicon
        #but check if there's actually a next replicon: is s < length(start)?
        #add coverage to Rlelist
        #switch to next replicon, repeat

        if (s < length(start)){
          chr_names <<- append(chr_names, df[start[s],]$V2)
          if (start[s] == 1 | start[s] == 2){

            if (is.null(list_cov)){
              rl4 <- IRanges(as.numeric(df[(start[s]+1):(start[s+1]-1), ]$V1),
                             as.numeric(df[(start[s]+1):(start[s+1]-1), ]$V1))

              scores <- as.numeric(df[(start[s]+1):(start[s+1]-1), ]$V2)

              #calculate coverage normally if the chunk contains a replicon
              cov <- abs(coverage(rl4, weight = scores))

              list_cov <<- cov

            }
            else{
              if(start[s] == 2){
                if (!(ch==0) & !(df[df$id==1,]$V1 == (coord+1))){
                  list_cov <<- c(list_cov,
                                 abs(coverage(IRanges(c((coord+1):(as.numeric(df[df$id==1,]$V1)-1)),c((coord+1):(as.numeric(df[df$id==1,]$V1)-1))),
                                              shift = -(coord), weight = rep(0, (as.numeric(df[df$id==1,]$V1) - coord)-1))))
                }
                rl4 <- IRanges(as.numeric(df[1:(start[s]-1), ]$V1),
                               as.numeric(df[1:(start[s]-1), ]$V1))

                scores <- as.numeric(df[1:(start[s]+1), ]$V2)

                cov <- abs(coverage(rl4, shift = -(as.numeric(df$V1[1])-1), weight = scores))
                list_cov <<- c(list_cov, cov)
              }


              #create Rlelist element
              if(isEmpty(list_rle)){
                list_rle <<- RleList(list_cov)
              }
              else{
                list_rle <<- append(list_rle, RleList(list_cov))
              }
              list_cov <<- NULL
              #create RLE for new replicon
              #chr_names <- append(chr_names, df[start[s],]$V2)
              rl4 <- IRanges(as.numeric(df[(start[s]+1):(start[s+1]-1), ]$V1),
                             as.numeric(df[(start[s]+1):(start[s+1]-1), ]$V1))

              scores <- as.numeric(df[(start[s]+1):(start[s+1]-1), ]$V2)
              cov <- abs(coverage(rl4, weight = scores))

              list_cov <<- cov
            }
          }

          else{#if start is further down in the chunk

            #for the first replicon
            #an RLE needs to be created from the beginning to update the prev replicon
            if (s==1){
              rl4 <- IRanges(as.numeric(df[1:(start[s]-1), ]$V1),
                             as.numeric(df[1:(start[s]-1), ]$V1))

              scores <- as.numeric(df[1:(start[s]-1), ]$V2)
              cov <- abs(coverage(rl4, shift = -(as.numeric(df$V1[1])-1), weight = scores))

              if (!(ch==0) & !(df[df$id==1,]$V1 == (coord+1))){
                cov <- c(abs(coverage(IRanges(c((coord+1):(as.numeric(df[df$id==1,]$V1)-1)),c((coord+1):(as.numeric(df[df$id==1,]$V1)-1))),
                                      shift = -(coord), weight = rep(0, (as.numeric(df[df$id==1,]$V1) - coord)-1))), cov)
              }


              list_cov <<- c(list_cov, cov)
            }

            #update RLElist anyway
            #if after first replicon, create and update rlelist
            #if not, just update rlelist and proceed to creating the new rle
            if(isEmpty(list_rle)){
              list_rle <<- RleList(list_cov)
            }
            else{
              list_rle <<- append(list_rle, RleList(list_cov))
            }
            list_cov <<- NULL

            #create RLE for new replicon
            if ((start[s]) < nrow(df)){
              rl4 <- IRanges(as.numeric(df[(start[s]+1):(start[s+1]-1), ]$V1),
                             as.numeric(df[(start[s]+1):(start[s+1]-1), ]$V1))

              scores <- as.numeric(df[(start[s]+1):(start[s+1]-1), ]$V2)
              cov <- abs(coverage(rl4, weight = scores))

              list_cov <<- cov
            }

          }
        }

        #if we've reached the final replicon in the chunk
        #RLE from rep position until the end
        else{
          #update RLElist
          if(isEmpty(list_rle)){
            list_rle <<- RleList(list_cov)
          }
          else{
            list_rle <<- append(list_rle, RleList(list_cov))
          }
          list_cov <<- NULL

          chr_names <<- append(chr_names, df[start[s],]$V2)

          if ((start[s]) < nrow(df)){
            rl4 <- IRanges(as.numeric(df[(start[s]+1):nrow(df), ]$V1),
                           as.numeric(df[(start[s]+1):nrow(df), ]$V1))

            scores <- as.numeric(df[(start[s]+1):nrow(df), ]$V2)
            cov <- abs(coverage(rl4, weight = scores))

            #assign(list_cov, cov, envir = .GlobalEnv)
            list_cov <<- cov
          }
        }
      }
    }
    #update coord to the last pos in the chunk
    #if it's a new rep, then variable/fixedStep is recorded
    #shouldn't matter...
    coord <<- as.numeric(df[df$id == nrow(df),]$V1)

  }

  #intialize
  list_cov <- NULL #list of coverages to be combined into an Rlelist object
  list_rle <- NULL #Rlelist object
  chr_names <- NULL #replicons
  coord <- 0 #stores the final position of the previous chunk to check if
  #the first position of the next chunk is exactly coord+1
  #else create rle of 0s
  f_wigs <- NULL
  n_rows <- 0
  chunks <- 0
  df_path <- getwd()
  print(df_path)

  wigs <- plyr::llply(data_table$file, function(f){

    #wigs.df <- csv_to_disk.frame(f, fill=T, header = F,
                                 #keep.rownames=T,
    #                             outdir = paste0(df_path,"/tmp_diskframe"))

    n_rows <<- nrow(fread(f, select = 1L, fill = T, header = F))
    chunks <- ceiling(n_rows/10000)


    list_cov <<- NULL
    list_rle <<- NULL
    chr_names <<- NULL
    coord <<- 0

    for (i in 0:(chunks-1)){
      if (i < (chunks-1)){

        if (i < 2){
          f_wigs[[i+1]] <- add_chunks(f,i, 1)
          covr <- calc_rle(f_wigs[[i+1]],i)

        }
        else if(i %% 2 == 0){
          even <- add_chunks(f,i, 1)
          covr <- calc_rle(even,i)#even
        }
        else if(i %% 2 == 1){
          odd <- add_chunks(f,i, 1)
          covr <- calc_rle(odd,i)#odd
        }


      }
      else{
        f_wigs[[i+1]] <- add_chunks(f,i, 2)
        covr <- calc_rle(f_wigs[[i+1]],i)
      }
    }
    #add the final RLE to the list
    list_rle <<- append(list_rle, RleList(list_cov))#or covr
    names(list_rle)<<-chr_names

    return(list_rle)

  }, .progress = "text")






  #####
  #check replicons identical, get lengths
  if(missing(replicons)){
    #cool, just comparing replicon_tags in the RLElist and
    #assigning names to replicons
    replicons <- names(wigs[[1]])
  }
  #print(replicons)


  #check if replicons extracted from expermient 1 are the same for all
  #the experiments
  plyr::l_ply(wigs, function(x){
    if(length(setdiff(replicons, names(x))) > 0){
      stop("Replicons differ between experiments")
    }
  })


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


  return(wigs)
}

