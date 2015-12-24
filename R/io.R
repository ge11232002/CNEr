### -----------------------------------------------------------------
### read the bed file (with only 3 columns) into GRanges.
### Exported!
readBed <- function(bedFile){
  ## GRanges: 1-based start
  ## bed file: 0-based start
  
  #bed <- .Call2("myReadBed", bedFile, PACKAGE="CNEr")
  bed <- read_tsv(bedFile, col_names=FALSE, comment="track")
  ## We only need the first three columns of the bed file, 
  ## but keep the strand information when available
  if(is.null(bed[[4]])){
    strands <- factor("+")
  }else{
    strands <- bed[[4]]
  }
  bed <- GRanges(seqnames=Rle(bed[[1]]),
                 ranges=IRanges(start=bed[[2]]+1L, end=bed[[3]]),
                 strand=strands)
  return(bed)
}

### -----------------------------------------------------------------
### read the axt files into an axt object.
### Exported!
readAxt <- function(axtFiles){
  # Read axt files into R axt object.
  # The coordinates are 1-based for start and end.
  index_noexists <- !file.exists(axtFiles)
  if(any(index_noexists)){
    stop("No such file ", paste(axtFiles[index_noexists], sep=" "))
  }
  myAxt <- .Call2("myReadAxt", axtFiles, PACKAGE="CNEr")
  axts <- Axt(targetRanges=GRanges(seqnames=Rle(myAxt[[1]]),
                                   ranges=IRanges(start=myAxt[[2]],
                                                  end=myAxt[[3]]),
                                   strand=Rle(myAxt[[4]])),
              targetSeqs=DNAStringSet(myAxt[[5]]),
              queryRanges=GRanges(seqnames=Rle(myAxt[[6]]),
                                  ranges=IRanges(start=myAxt[[7]],
                                                 end=myAxt[[8]]),
                                  strand=Rle(myAxt[[9]])),
              querySeqs=DNAStringSet(myAxt[[10]]),
              score=myAxt[[11]],
              symCount=myAxt[[12]]
            )
  return(axts)
}

### -----------------------------------------------------------------
### read the axt files and return the widths of all the alignments
### Exported!
axtInfo <- function(axtFiles){
  ans <- .Call2("axt_info", axtFiles, PACKAGE="CNEr")
  return(ans)
}

### -----------------------------------------------------------------
### write the Axt object to an axt file
### Exported!
writeAxt <- function(axt, con){
  firstLine <- paste(0:(length(axt)-1), seqnames(targetRanges(axt)),
                     start(targetRanges(axt)), end(targetRanges(axt)),
                     seqnames(queryRanges(axt)),
                     start(queryRanges(axt)), end(queryRanges(axt)),
                     strand(queryRanges(axt)), score(axt)
                     )
  secondLine <- targetSeqs(axt)
  thirdLine <- querySeqs(axt)
  wholeLines <- paste(firstLine, as.character(targetSeqs(axt)), 
                      as.character(querySeqs(axt)),
                      "", sep="\n")
  writeLines(wholeLines, con)
}



