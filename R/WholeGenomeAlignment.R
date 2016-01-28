### -----------------------------------------------------------------
### scoringMatrix
### Exported!
scoringMatrix <- function(distance=c("far", "medium", "close")){
  distance <- match.arg(distance)
  lastzMatrix <- list(medium=matrix(c(91, -114, -31, -123,
                                      -114, 100, -125, -31,
                                      -31, -125, 100,-114,
                                      -123, -31, -114, 91),
                                    nrow=4, ncol=4,
                                    dimnames=list(c("A", "C", "G", "T"),
                                                  c("A", "C", "G", "T"))
  ),
  far=matrix(c(91, -90, -25, -100,
               -90, 100, -100, -25,
               -25, -100, 100, -90,
               -100, -25, -90, 91),
             nrow=4, ncol=4,
             dimnames=list(c("A", "C", "G", "T"),
                           c("A", "C", "G", "T"))
  ),
  near=matrix(c(90, -330, -236, -356,
                -330, 100, -318, -236,
                -236, -318, 100, -330,
                -356, -236, -330, 90),
              nrow=4, ncol=4,
              dimnames=list(c("A", "C", "G", "T"),
                            c("A", "C", "G", "T"))
  )
  )
  return(lastzMatrix[[distance]])
}


### -----------------------------------------------------------------
### lastz wrapper
### Exported!
lastz <- function(assemblyTarget, assemblyQuery,
                  outputDir=".",
                  chrsTarget=NULL, chrsQuery=NULL, 
                  distance=c("far", "medium", "close"),
                  binary="lastz",
                  mc.cores=getOption("mc.cores", 2L),
                  echoCommand=FALSE){
  distance <- match.arg(distance)
  if(!all(file.exists(c(assemblyTarget, assemblyQuery)))){
    stop(assemblyTarget, " and ", assemblyQuery, " must exist!")
  }
  if(file_ext(assemblyTarget) != "2bit" || file_ext(assemblyQuery) != "2bit"){
    stop("The assembly must be in .2bit format!")
  }
  
  matrixFile <- tempfile(fileext=".lastzMatrix")
  ## The options used here is taken from RunLastzChain_sh.txt 
  ## genomewiki.ucsc.edu.
  ## http://genomewiki.ucsc.edu/images/9/93/RunLastzChain_sh.txt.
  
  ## B=0, --strand=plus: Search the forward strand only 
  ## (the one corresponding to the query specifier).
  ## By default, both strands are searched.
  
  ## C=0, --nochain: Skip the chaining stage. 
  ## By default, the chaining stage is skipped.
  
  ## E=30,150;O=600,400: Set the score penalties for opening and 
  ## extending a gap.
  
  ## H=0, 2000; --inner=<score>: Perform additional alignment
  ## between the gapped alignment blocks,
  ## using (presumably) more sensitive alignment parameters.
  ## By default this is not performed.
  
  ## K=4500,3000,2200; --hspthresh=<score>:
  ## Set the score threshold for the x-drop extension method;
  ## HSPs scoring lower are discarded. By default, use the entropy adjustment.
  
  ## L=3000,6000; --gappedthresh=<score>: 
  ## Set the threshold for gapped extension;
  ## alignments scoring lower than <score> are discarded.
  ## By default gapped extension is performed,
  ## and alignment ends are trimmed to the locations giving the maximum score.
  
  ## M=254,50,50;--masking=<count>:  
  ## Dynamically mask the target sequence by 
  ## excluding any positions that appear in too many alignments from 
  ## further consideration for seeds. 
  ## By default, a step of 1 is used, 
  ## no words are removed from the target seed word position table,
  ## dynamic masking is not performed,
  ## and no target capsule or segment file is used.
  
  ## T=1,2;--seed=12of19:  Seeds require a 19-bp word 
  ## with matches in 12 specific positions (1110100110010101111).
  ## By default the 12-of-19 seed is used,
  ## one transition is allowed (except with quantum DNA),
  ## the hits are not filtered, twins are not required,
  ## and hash collisions are not recovered.
  
  ##  Y=15000,9400,3400; --ydrop=<dropoff>:
  ## Set the threshold for terminating gapped extension;
  ## this restricts the endpoints of each local alignment
  ## by limiting the local region around each anchor 
  ## in which extension is performed.
  
  lastzOptions <- list(
    near=paste0("C=0 E=150 H=0 K=4500 L=3000 M=254 O=600 T=2 Y=15000 Q=",
                matrixFile),
    medium=paste0("C=0 E=30 H=0 K=3000 L=3000 M=50 O=400 T=1 Y=9400",
                  matrixFile),
    far=paste0("C=0 E=30 H=2000 K=2200 L=6000 M=50 O=400 T=2 Y=3400 Q=", matrixFile)
  )
  write.table(scoringMatrix(distance), file=matrixFile, quote=FALSE,
              sep=" ", row.names=FALSE, col.names=TRUE)
  message("Starting lastz")
  if(is.null(chrsTarget)){
    chrsTarget <- seqnames(seqinfo(TwoBitFile(assemblyTarget)))
  }else{
    stopifnot(all(chrsTarget %in% 
                    seqnames(seqinfo(TwoBitFile(assemblyTarget)))))
  }
  if(is.null(chrsQuery)){
    chrsQuery <- seqnames(seqinfo(TwoBitFile(assemblyQuery)))
  }else{
    stopifnot(all(chrsQuery %in% seqnames(seqinfo(TwoBitFile(assemblyQuery)))))
  }
  outputToReturn <- c()
  format <- "lav"
  .runLastz <- function(chrTarget, chrQuery, assemblyTarget, assemblyQuery,
                        format="lav"){
    ## Deal the "|" in chr name
    output <- file.path(outputDir,
                        paste0(gsub("|", "\\|", chrTarget, fixed=TRUE),
                        ".", sub("\\..*$", "", basename(assemblyTarget)),
                        "-", gsub("|", "\\|", chrQuery, fixed=TRUE),
                        ".", sub("\\..*$", "", basename(assemblyQuery)),
                        ".", format))
    cmd <- paste0(binary, " ", assemblyTarget, "/",
                  gsub("|", "\\|", chrTarget, fixed=TRUE), " ",
                  assemblyQuery, "/",
                  gsub("|", "\\|", chrQuery, fixed=TRUE), " ",
                  lastzOptions[[distance]],
                  " --format=", format,
                  " --output=", output,
                  " --markend")
    if(echoCommand){
      output <- cmd
    }else{
      if(!file.exists(output)){
        my.system(cmd)
      }else{
        warning("The output ", output, " already exists! Skipping..")
      }
    }
    return(output)
  }
  combinations <- expand.grid(chrsTarget, chrsQuery,
                              stringsAsFactors=FALSE)
  mc.cores <- min(mc.cores, nrow(combinations))
  ans <- mcmapply(.runLastz, combinations[[1]], combinations[[2]],
                  MoreArgs=list(assemblyTarget=assemblyTarget,
                                assemblyQuery=assemblyQuery,
                                format=format),
                  mc.cores=mc.cores)
  unlink(matrixFile)
  invisible(unname(ans))
}

### -----------------------------------------------------------------
### validateLastz: validate the lav file generated from lastz is complete.
### Only works on unix-like system. tail is the efficient way.
### Not Exported!
validateLastz <- function(lavs){
  filesNotCompleted <- c()
  for(lav in lavs){
    lastLine <- my.system(paste("tail -n 1", lav), intern=TRUE)
    if(lastLine != "# lastz end-of-file"){
      filesNotCompleted <- c(filesNotCompleted, lav)
    }
  }
  return(filesNotCompleted)
}

### -----------------------------------------------------------------
### lavToPsl: convert lav files to psl files
### Exported!
lavToPsl <- function(lavs, 
                     psls=sub("\\.lav$", ".psl", lavs, ignore.case=TRUE),
                     removeLav=TRUE){
  .runLav2Psl <- function(lav, psl){
    arguments <- c(lav, psl)
    system2(command="lavToPsl", args=arguments)
    return(psl)
  }
  ans <- mapply(.runLav2Psl, lavs, psls)
  if(removeLav){
    unlink(lavs)
  }
  invisible(unname(ans))
}

### -----------------------------------------------------------------
### axtChain: wrapper for axtChain
### Exported!
axtChain <- function(psls,
                     chains=sub("\\.psl$", ".chain", psls, ignore.case=TRUE), 
                     assemblyTarget, assemblyQuery,
                     distance=c("far", "medium", "far"),
                     removePsl=TRUE){
  distance <- match.arg(distance)
  
  if(file_ext(assemblyTarget) != "2bit" || file_ext(assemblyQuery) != "2bit"){
    stop("The assembly must be in .2bit format!")
  }
  
  chainOptions <- list(near="-minScore=5000 -linearGap=medium",
                       medium="-minScore=3000 -linearGap=medium",
                       far="-minScore=5000 -linearGap=loose")
  matrixFile <- tempfile(fileext=".lastzMatrix")
  write.table(scoringMatrix(distance), file=matrixFile, quote=FALSE,
              sep=" ", row.names=FALSE, col.names=TRUE)
  if(distance == "far"){
    linearGap <- "loose"
  }else{
    linearGap <- "medium"
  }
  
  .runAxtChain <- function(psl, chain){
    arguments <- c("-psl", chainOptions[[distance]],
                   paste0("-scoreScheme=", matrixFile),
                   paste0("-linearGap", linearGap),
                   psl, assemblyTarget, assemblyQuery, chain)
    system2(command="axtChain", args=arguments)
    return(chain)
  }
  ans <- mapply(.runAxtChain, psls, chains)
  
  if(removePsl){
    unlink(psls)
  }
  
  unlink(matrixFile)
  invisible(unname(ans))
}
