### -----------------------------------------------------------------
### read the bed file into GRanges.
### Exported!
readBed <- function(bedFile){
  ## GRanges: 1-based start
  ## bed file: 0-based start
  
  #bed <- .Call2("myReadBed", bedFile, PACKAGE="CNEr")
  bed <- read_tsv(bedFile, col_names=FALSE, comment="track")
  ## We only need the first three columns of the bed file, 
  ## but keep the strand information when available
  if(ncol(bed) == 3L){
    strands <- factor("+")
  }else{
    strands <- bed[[6]]
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
  axts <- Axt(targetRanges=GRanges(seqnames=myAxt[[1]],
                                   ranges=IRanges(start=myAxt[[2]],
                                                  end=myAxt[[3]]),
                                   strand=myAxt[[4]]),
              targetSeqs=DNAStringSet(myAxt[[5]]),
              queryRanges=GRanges(seqnames=myAxt[[6]],
                                  ranges=IRanges(start=myAxt[[7]],
                                                 end=myAxt[[8]]),
                                  strand=myAxt[[9]]),
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
  index_noexists <- !file.exists(axtFiles)
  if(any(index_noexists)){
    stop("No such file ", paste(axtFiles[index_noexists], sep=" "))
  }
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

### -----------------------------------------------------------------
### read RepeatMasker out file into a GRanges object
### Exported!
read.rmMask.GRanges <- function(fn){
  rmMaskOut <- read.table(fn, header=FALSE, sep="", skip=3, as.is=TRUE,
                          col.names=1:16, fill=TRUE)
  rmMaskGRanges <- GRanges(seqnames=rmMaskOut$X5,
                           ranges=IRanges(start=rmMaskOut$X6,
                                          end=rmMaskOut$X7),
                           strand=ifelse(rmMaskOut$X9=="+", "+", "-"),
                           name=rmMaskOut$X10,
                           type=rmMaskOut$X11,
                           score=rmMaskOut$X1)
  return(rmMaskGRanges)
}

### -----------------------------------------------------------------
### save the CNE tables into a local SQLite database
### Exported!!
saveCNEToSQLite <- function(x, dbName, tableName=NULL, overwrite=FALSE){
  ## by default tableName is in the format "danRer7_hg19_49_50"
  if(is.null(tableName)){
    tableName <- paste(sub("\\.2bit", "", basename(x@assembly1Fn)),
                       sub("\\.2bit", "", basename(x@assembly2Fn)),
                       x@identity, x@window, sep="_")
  }
  if(length(CNEFinal(x)) == 0L){
    warning("There is no CNEs.")
  }
  cneFinal <- as.data.frame(CNEFinal(x))
  colnames(cneFinal) <- sub("^X\\.", "", colnames(cneFinal))
  
  ## Create the bin column
  cneFinal$first.bin <- binFromCoordRange(cneFinal$first.start,
                                          cneFinal$first.end)
  cneFinal$last.bin <- binFromCoordRange(cneFinal$last.start,
                                         cneFinal$last.end)
  cneFinal <- cneFinal[ ,c("first.bin", "first.seqnames",
                           "first.start", "first.end",
                           "last.bin",
                           "last.seqnames", "last.start", "last.end"
  )]
  
  ## SQLite
  con <- dbConnect(SQLite(), dbname=dbName)
  on.exit(dbDisconnect(con))
  dbWriteTable(con, tableName, cneFinal, row.names=FALSE,
               overwrite=overwrite)
  invisible(tableName)
}