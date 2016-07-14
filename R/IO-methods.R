### -----------------------------------------------------------------
### seqinfoFn: get the Seqinfo object from fasta or twoBit file.
### Not Exported!
seqinfoFn <- function(fn){
  fileType <- file_ext(fn)
  genome <- sub("\\..*$", "", basename(fn))
  if(fileType %in% c("fa", "fasta")){
    # fasta file
    seqlengths1 <- fasta.seqlengths(fn)
    names(seqlengths1) <- sapply(strsplit(names(seqlengths1), " "), "[", 1)
    ## Only keep the first part of fasta lines.
    ans <- Seqinfo(seqnames=names(seqlengths1),
                   seqlengths=seqlengths1,
                   genome=genome)
  }else if(fileType == "2bit"){
    # 2bit file
    ans <- seqinfo(TwoBitFile(fn))
  }
  genome(ans) <- genome
  return(ans)
}

### -----------------------------------------------------------------
### readBed: read the bed file into GRanges.
### Exported!
readBed <- function(bedFile, assemblyFn=NULL){
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
  seqinfoBed <- NULL
  if(!is.null(assemblyFn)){
    seqinfoBed <- seqinfoFn(assemblyFn)
  }
  ans <- GRanges(seqnames=Rle(bed[[1]]),
                 ranges=IRanges(start=bed[[2]]+1L, end=bed[[3]]),
                 strand=strands, seqinfo=seqinfoBed)
  return(ans)
}

### -----------------------------------------------------------------
### read the axt files into an axt object.
### Exported!
readAxt <- function(axtFiles, tAssemblyFn=NULL, qAssemblyFn=NULL){
  # Read axt files into R axt object.
  # The coordinates are 1-based for start and end.
  index_noexists <- !file.exists(axtFiles)
  if(any(index_noexists)){
    stop("No such file ", paste(axtFiles[index_noexists], sep=" "))
  }
  
  # Prepare the seqinfo when available
  seqinfoTarget <- NULL
  if(!is.null(tAssemblyFn)){
    seqinfoTarget <- seqinfoFn(tAssemblyFn)
  }
  seqinfoQuery <- NULL
  if(!is.null(qAssemblyFn)){
    seqinfoQuery <- seqinfoFn(qAssemblyFn)
  }
  
  ## Extend the absolute paths of files
  axtFiles <- normalizePath(axtFiles)
  myAxt <- .Call2("myReadAxt", axtFiles, PACKAGE="CNEr")
  axts <- Axt(targetRanges=GRanges(seqnames=myAxt[[1]],
                                   ranges=IRanges(start=myAxt[[2]],
                                                  end=myAxt[[3]]),
                                   strand=myAxt[[4]],
                                   seqinfo=seqinfoTarget),
              targetSeqs=DNAStringSet(myAxt[[5]]),
              queryRanges=GRanges(seqnames=myAxt[[6]],
                                  ranges=IRanges(start=myAxt[[7]],
                                                 end=myAxt[[8]]),
                                  strand=myAxt[[9]],
                                  seqinfo=seqinfoQuery),
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
### save the CNE class or GRangePairs object into a local SQLite database
### Exported!!
saveCNEToSQLite <- function(x, dbName, tableName=NULL, overwrite=FALSE){
  ## by default tableName is in the format "danRer7_hg19_49_50"
  if(is.null(tableName)){
    tableName <- paste(sub("\\.2bit", "", basename(x@assembly1Fn)),
                       sub("\\.2bit", "", basename(x@assembly2Fn)),
                       x@identity, x@window, sep="_")
  }
  if(class(x) == "CNE"){
    ## CNE class
    cneFinal <- as.data.frame(CNEFinal(x))
  }else if(class(x) == "GRangePairs"){
    ## GRangePairs class 
    cneFinal <- as.data.frame(x)
  }else{
    stop(" x must be a CNE class or GRangePairs class.")
  }
  if(nrow(cneFinal) == 0L){
    warning("There is no CNEs.")
  }
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

### -----------------------------------------------------------------
### read CNE from a local SQLite database
### Exported!
readCNERangesFromSQLite <- function(dbName, tableName,
                                    chr=NULL, start=NULL, end=NULL,
                                    whichAssembly=c("first", "last"),
                                    minLength=NULL){
  nrGraphs <- 1
  ## Let's make nrGraphs=1, make all the cnes together.
  whichAssembly <- match.arg(whichAssembly)
  con <- dbConnect(SQLite(), dbname=dbName)
  on.exit(dbDisconnect(con))
  
  if(is.null(chr) && is.null(start) & is.null(end)){
    # 1. fetch all the CNEs: chr=NULL, start=NULL, end=NULl
    #sqlCmd <- switch(whichAssembly,
    #  "first"=paste("SELECT [first.seqnames],[first.start],[first.end] from",
    #                tableName),
    #  "last"=paste("SELECT [last.seqnames],[last.start],[last.end] from",
    #               tableName),
    #  "all"=paste("SELECT [first.seqnames],[first.start],[first.end]",
    #              "[last.seqnames],[last.start],[last.end] from",
    #              tableName)
    #  )
    sqlCmd <- paste("SELECT [first.seqnames],[first.start],[first.end],[last.seqnames],[last.start],[last.end] from",
                    tableName)
  }else if(!is.null(chr) && is.null(start) && is.null(end)){
    # 2. fetch all CNEs on chromosomes chr
    sqlCmd <- switch(whichAssembly,
      "first"=paste("SELECT [first.seqnames],[first.start],[first.end],[last.seqnames],[last.start],[last.end] from",
                    tableName, "WHERE [first.seqnames] IN (",
                    paste(paste0("'", chr, "'"), collapse=","),
                    ")"),
      "last"=paste("SELECT [first.seqnames],[first.start],[first.end],[last.seqnames],[last.start],[last.end] from",
                   tableName, "WHERE [last.seqnames] IN (",
                   paste(paste0("'", chr, "'"), collapse=","),
                   ")")
    )
  }else if(!is.null(chr) && !is.null(start) && !is.null(end)){
    # 3. fetch all CNEs on potentially multiple chr, start, end
    CNEstart <- as.integer(start)
    CNEend <- as.integer(end)
    sqlCmd <- switch(whichAssembly,
      "first"=paste("SELECT [first.seqnames],[first.start],[first.end],[last.seqnames],[last.start],[last.end] from",
                    tableName, "WHERE", paste("([first.seqnames]=",
                                              paste0("'", chr, "'"),
                                              "AND [first.start] >=", CNEstart,
                                              "AND [first.end] <=",
                                              CNEend, "AND",
                      binRestrictionString(CNEstart, CNEend,
                                           "[first.bin]"),
                      ")", collapse=" OR ")),
      "last"=paste("SELECT [first.seqnames],[first.start],[first.end],[last.seqnames],[last.start],[last.end] from",
                   tableName, "WHERE", paste("([last.seqnames]=",
                                             paste0("'", chr, "'"),
                                             "AND [last.start] >=", CNEstart,
                                             "AND [last.end] <=",
                                             CNEend, "AND",
                     binRestrictionString(CNEstart, CNEend,
                                          "[last.bin]"),
                     ")", collapse=" OR "))
      )
  }else{
    stop("Unsupported search criteria!")
  }
  if(!is.null(minLength))
    sqlCmd <- paste(sqlCmd, "AND [first.end]-[first.start]+1 >=", minLength, 
                    "AND [last.end]-[last.start]+1 >=", minLength)
  fetchedCNE <- dbGetQuery(con, sqlCmd)
  firstGRanges <- GRanges(seqnames=fetchedCNE[ ,1],
                          ranges=IRanges(start=fetchedCNE[ ,2],
                                         end=fetchedCNE[,3]),
                          strand="*")
  lastGRanges <- GRanges(seqnames=fetchedCNE[ ,4],
                         ranges=IRanges(start=fetchedCNE[ ,5],
                                        end=fetchedCNE[ ,6]),
                         strand="*")
  ans <- GRangePairs(first=firstGRanges, last=lastGRanges)
  # if(nrGraphs == 1){
  #   sqlCmd <- switch(whichAssembly,
  #                    "first"=paste("SELECT [first.start],[first.end] from",
  #                                  tableName, 
  #                                  "WHERE [first.seqnames]=", 
  #                                  paste0("'", chr, "'"), 
  #                                  "AND [first.start] >=", CNEstart, 
  #                                  "AND [first.end] <=", 
  #                                  CNEend, "AND", 
  #                                  binRestrictionString(CNEstart, CNEend,
  #                                                       "[first.bin]")),
  #                    "last"=paste("SELECT [last.start],[last.end] from",
  #                                 tableName, 
  #                                 "WHERE [last.seqnames]=", 
  #                                 paste0("'", chr, "'"), 
  #                                 "AND [last.start] >=", CNEstart, 
  #                                 "AND [last.end] <=", 
  #                                 CNEend, "AND", 
  #                                 binRestrictionString(CNEstart, CNEend, 
  #                                                      "[last.bin]"))
  #   )
    
  # }else if(nrGraphs > 1){
  #   ## This chunk is not executed for now.
  #   sqlCmd <- switch(whichAssembly,
  #                    "L"=paste("SELECT chr2,start1,end1 from", tableName, 
  #                              "WHERE chr1=", paste0("'", chr, "'"), 
  #                              "AND start1 >=", CNEstart, "AND end1 <=", 
  #                              CNEend, "AND", 
  #                              binRestrictionString(CNEstart, CNEend, "bin1")),
  #                    "R"=paste("SELECT chr1,start2,end2 from", tableName, 
  #                              "WHERE chr2=", paste0("'", chr, "'"), 
  #                              "AND start2 >=", CNEstart, "AND end2 <=", 
  #                              CNEend, "AND", 
  #                              binRestrictionString(CNEstart, CNEend, "bin2"))
  #   )
  #   if(!is.null(minLength))
  #     sqlCmd <- paste(sqlCmd, "AND end1-start1+1 >=", minLength, 
  #                     "AND end2-start2+1 >=", minLength)
  #   fetchedCNE <- dbGetQuery(con, sqlCmd)
  #   fetchedCNE <- GRanges(seqnames=fetchedCNE[ ,1], 
  #                         ranges=IRanges(start=fetchedCNE[ ,2], 
  #                                        end=fetchedCNE[ ,3]))
  # }
  return(ans)
}