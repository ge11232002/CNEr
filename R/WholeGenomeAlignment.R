### -----------------------------------------------------------------
### scoringMatrix
### Exported!
scoringMatrix <- function(distance=c("far", "medium", "near")){
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
                  distance=c("far", "medium", "near"),
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
  if(!dir.exists(outputDir)){
    dir.create(outputDir, recursive=TRUE)
  }
  
  matrixFile <- tempfile(fileext=".lastzMatrix")
  on.exit(unlink(matrixFile))
  write.table(scoringMatrix(distance), file=matrixFile, quote=FALSE,
              sep=" ", row.names=FALSE, col.names=TRUE)
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
                     removeLav=TRUE, binary="lavToPsl"){
  .runLav2Psl <- function(lav, psl){
    arguments <- c(lav, psl)
    system2(command=binary, args=arguments)
    return(psl)
  }
  message("Run lavToPsl...")
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
                     distance=c("far", "medium", "near"),
                     removePsl=TRUE,
                     binary="axtChain"){
  distance <- match.arg(distance)
  
  if(file_ext(assemblyTarget) != "2bit" || file_ext(assemblyQuery) != "2bit"){
    stop("The assembly must be in .2bit format!")
  }
  
  chainOptions <- list(near="-minScore=5000 -linearGap=medium",
                       medium="-minScore=3000 -linearGap=medium",
                       far="-minScore=5000 -linearGap=loose")
  matrixFile <- tempfile(fileext=".lastzMatrix")
  on.exit(unlink(matrixFile))
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
    system2(command=binary, args=arguments)
    return(chain)
  }
  ans <- mapply(.runAxtChain, psls, chains)
  
  if(removePsl){
    unlink(psls)
  }
  
  invisible(unname(ans))
}

### -----------------------------------------------------------------
### chainMergeSort: Combine sorted files into larger sorted file.
### Exported!
chainMergeSort <- function(chains, assemblyTarget, assemblyQuery,
                           allChain=paste0(sub("\\.2bit$", "",
                                              basename(assemblyTarget),
                                              ignore.case=TRUE),
                                          ".",
                                          sub("\\.2bit$", "",
                                              basename(assemblyQuery),
                                              ignore.case=TRUE),
                                          ".all.chain"),
                           removeChains=TRUE,
                           binary="chainMergeSort"){
  chainFile <- tempfile(fileext=".chainMergeSort")
  on.exit(file.remove(chainFile))
  writeLines(chains, chainFile)
  
  arguments <- c(paste0("-inputList=", chainFile), paste(">", allChain))
  message("Run chainMergeSort...")
  system2(command=binary, args=arguments)
  
  if(removeChains){
    unlink(chains)
  }
  invisible(allChain)
}


### -----------------------------------------------------------------
### chainPreNet: Remove chains that don't have a chance of being netted
### Exported!
chainPreNet <- function(allChain, assemblyTarget, assemblyQuery,
                        allPreChain=paste0(sub("\\.2bit$", "",
                                               basename(assemblyTarget),
                                               ignore.case=TRUE),
                                          ".",
                                          sub("\\.2bit$", "",
                                              basename(assemblyQuery),
                                              ignore.case=TRUE),
                                          ".all.pre.chain"),
                        removeAllChain=TRUE, binary="chainPreNet"){
  
  ## prepare the genome size file
  target.sizesFile <- tempfile(fileext=".target.sizes")
  write.table(as.data.frame(seqinfo(
                TwoBitFile(assemblyTarget)))[ ,c("seqlengths"), drop=FALSE],
              file=target.sizesFile, row.names=TRUE, col.names=FALSE,
              quote=FALSE)
  query.sizesFile <- tempfile(fileext=".query.sizes")
  write.table(as.data.frame(seqinfo(
    TwoBitFile(assemblyQuery)))[ ,c("seqlengths"), drop=FALSE],
    file=query.sizesFile, row.names=TRUE, col.names=FALSE,
    quote=FALSE)
  on.exit(unlink(c(target.sizesFile, query.sizesFile)))
  
  ## run chainPreNet
  arguments <- c(allChain, target.sizesFile, query.sizesFile, allPreChain)
  message("Run chainPreNet...")
  system2(command=binary, args=arguments)
  
  if(removeAllChain){
    unlink(allChain)
  }
  invisible(allPreChain)
}

### -----------------------------------------------------------------
### chainNetSyntenic: Make alignment nets out of chains and 
### Add synteny info to net.
### Exported!
chainNetSyntenic <- function(allPreChain, assemblyTarget, assemblyQuery,
                             netSyntenicFile=paste0(
                               sub("\\.2bit$", "", basename(assemblyTarget),
                                   ignore.case=TRUE), ".",
                               sub("\\.2bit$", "", basename(assemblyQuery),
                                   ignore.case=TRUE), ".noClass.net"),
                             binaryChainNet="chainNet",
                             binaryNetSyntenic="netSyntenic"){
  ## prepare the genome size file
  target.sizesFile <- tempfile(fileext=".target.sizes")
  write.table(as.data.frame(seqinfo(
    TwoBitFile(assemblyTarget)))[ ,c("seqlengths"), drop=FALSE],
    file=target.sizesFile, row.names=TRUE, col.names=FALSE,
    quote=FALSE)
  query.sizesFile <- tempfile(fileext=".query.sizes")
  write.table(as.data.frame(seqinfo(
    TwoBitFile(assemblyQuery)))[ ,c("seqlengths"), drop=FALSE],
    file=query.sizesFile, row.names=TRUE, col.names=FALSE,
    quote=FALSE)
  on.exit(unlink(c(target.sizesFile, query.sizesFile)))
  
  ## chainNet
  message("Run chainNet...")
  target.net <- tempfile(fileext=".target.net")
  query.net  <- tempfile(fileext=".query.net")
  on.exit(unlink(c(target.net, query.net)))
  arguments <- c(allPreChain, target.sizesFile, query.sizesFile,
                 target.net, query.net)
  system2(command=binaryChainNet, args=arguments)
  
  ## netSyntenic
  message("Run netSyntenic...")
  arguments <- c(target.net, netSyntenicFile)
  system2(command=binaryNetSyntenic, args=arguments)
  
  invisible(netSyntenicFile)
}

### -----------------------------------------------------------------
### netToAxt: Convert net (and chain) to axt, and sort axt files
### Exported!
netToAxt <- function(in.net, in.chain, assemblyTarget, assemblyQuery,
                     axtFile=paste0(sub("\\.2bit$", "",
                                        basename(assemblyTarget),
                                        ignore.case=TRUE),
                                    ".",
                                    sub("\\.2bit$", "",
                                        basename(assemblyQuery),
                                        ignore.case=TRUE),
                                    ".net.axt"),
                     removeFiles=FALSE,
                     binaryNetToAxt="netToAxt", binaryAxtSort="axtSort"){
  
  ## netToAxt
  unsortedAxt <- tempfile(fileext=".unsorted.Axt")
  on.exit(file.remove(unsortedAxt))
  arguments <- c(in.net, in.chain, assemblyTarget, assemblyQuery, unsortedAxt)
  message("Run netAxt...")
  system2(command=binaryNetToAxt, args=arguments)
  
  ## axtSort
  arguments <- c(unsortedAxt, axtFile)
  message("Run axtSort...")
  system2(command=binaryAxtSort, args=arguments)
  
  ## Clean
  if(removeFiles){
    unlink(c(in.net, in.chain))
  }
  
  invisible(axtFile)
}

### -----------------------------------------------------------------
### last: wrapper function of last aligner
### Exported!
lastal <- function(db, queryFn,
                   outputFn=sub("\\.(fa|fasta)$", ".maf", 
                                paste(basename(db), basename(queryFn), sep=","),
                                ignore.case=TRUE),
                   distance=c("far", "medium", "near"),
                   binary="lastal",
                   mc.cores=getOption("mc.cores", 2L),
                   echoCommand=FALSE){
  distance <- match.arg(distance)
  
  if(file_ext(outputFn) != "maf" ){
    stop("The outputFn must be in .maf format!")
  }
  
  matrixFile <- tempfile(fileext=".lastzMatrix")
  on.exit(unlink(matrixFile))
  write.table(scoringMatrix(distance), file=matrixFile, quote=FALSE,
              sep=" ", row.names=TRUE, col.names=TRUE)
  
  ## -a: Gap existence cost.
  ## -b: Gap extension cost.
  ## -e: Minimum alignment score.
  ## -p: Specify a match/mismatch score matrix.  Options -r and -q will be ignored.
  ## -s: Specify which query strand should be used: 0 means reverse only, 1 means forward only, and 2 means both.
  lastOptiosn <- list(near=paste("-a 600 -b 150 -e 3000 -p",
                                 matrixFile, "-s 2"),
                      medium=paste("-a 400 -b 30 -e 4500 -p",
                                   matrixFile, "-s 2"),
                      far=paste("-a 400 -b 30 -e 6000 -p",
                                matrixFile, "-s 2")
                      )
  if(mc.cores == 1L){
    cmd <- paste("lastal", lastOptiosn[[distance]],
                 "-f 1",
                 db, queryFn, ">", outputFn)
  }else{
    message("lastal require `parallel` installed on the machine to run in parallel,
            otherwise it will fail!")
    cmd <- paste("parallel-fasta", "-j", mc.cores, "--compress",
                 "\"lastal", lastOptiosn[[distance]],
                 "-f 1", db, "\"", "<", queryFn,
                 ">", outputFn)
  }
  
  message("Run lastal..")
  if(echoCommand){
    message(cmd)
  }else{
    my.system(cmd)
  }
  invisible(outputFn)
}