
### -----------------------------------------------------------------
### The main function for scanning the axts and get the CNEs
### Not exported!
ceScanR <- function(axts, tFilter=NULL, qFilter=NULL, qSizes=NULL, 
                   thresholds=c("49_50")){
  ## Here the returned tStart and qStart are 1-based coordinates. 
  ## Of course ends are also 1-based.
  if(!is.null(qFilter))
    if(is.null(qSizes) || !is(qSizes, "Seqinfo"))
      stop("qSizes must exist and be a Seqinfo object when qFilter exists")
  
  winSize <- as.integer(sapply(strsplit(thresholds, "_"), "[", 2))
  minScore <- as.integer(sapply(strsplit(thresholds, "_"), "[", 1))
  resFiles <- tempfile(pattern=paste(minScore, winSize, "ceScan", sep="-"), 
                       tmpdir=tempdir(), fileext="")
  ## How stupid I am...
  if(is.null(tFilter) && is.null(qFilter)){
    .Call2("myCeScan", NULL, NULL, NULL,
           NULL, NULL, NULL,
           NULL, NULL,
           as.character(seqnames(targetRanges(axts))),
           start(targetRanges(axts)), end(targetRanges(axts)),
           as.character(strand(targetRanges(axts))),
           as.character(targetSeqs(axts)),
           as.character(seqnames(queryRanges(axts))),
           start(queryRanges(axts)), end(queryRanges(axts)),
           as.character(strand(queryRanges(axts))),
           as.character(querySeqs(axts)),
           score(axts), symCount(axts), winSize, minScore,
           as.character(resFiles),
           PACKAGE="CNEr")
  }else if(is.null(tFilter) && !is.null(qFilter)){
    .Call2("myCeScan", NULL, NULL, NULL,
           as.character(seqnames(qFilter)), start(qFilter), end(qFilter),
           as.character(seqnames(qSizes)), as.integer(seqlengths(qSizes)),
           as.character(seqnames(targetRanges(axts))),
           start(targetRanges(axts)), end(targetRanges(axts)),
           as.character(strand(targetRanges(axts))),
           as.character(targetSeqs(axts)),
           as.character(seqnames(queryRanges(axts))),
           start(queryRanges(axts)), end(queryRanges(axts)),
           as.character(strand(queryRanges(axts))),
           as.character(querySeqs(axts)),
           score(axts), symCount(axts), winSize, minScore,
           as.character(resFiles),
           PACKAGE="CNEr")
  }else if(!is.null(tFilter) && is.null(qFilter)){
    .Call2("myCeScan", as.character(seqnames(tFilter)), 
           start(tFilter), end(tFilter),
           NULL, NULL, NULL,
           NULL, NULL,
           as.character(seqnames(targetRanges(axts))),
           start(targetRanges(axts)), end(targetRanges(axts)),
           as.character(strand(targetRanges(axts))),
           as.character(targetSeqs(axts)),
           as.character(seqnames(queryRanges(axts))),
           start(queryRanges(axts)), end(queryRanges(axts)),
           as.character(strand(queryRanges(axts))),
           as.character(querySeqs(axts)),
           score(axts), symCount(axts), winSize, minScore,
           as.character(resFiles),
           PACKAGE="CNEr")
  }else{
    .Call2("myCeScan", as.character(seqnames(tFilter)), 
           start(tFilter), end(tFilter),
              as.character(seqnames(qFilter)), start(qFilter), end(qFilter),
              as.character(seqnames(qSizes)), as.integer(seqlengths(qSizes)), 
              as.character(seqnames(targetRanges(axts))), 
              start(targetRanges(axts)), end(targetRanges(axts)), 
              as.character(strand(targetRanges(axts))), 
              as.character(targetSeqs(axts)),
              as.character(seqnames(queryRanges(axts))), 
              start(queryRanges(axts)), end(queryRanges(axts)), 
              as.character(strand(queryRanges(axts))), 
              as.character(querySeqs(axts)),
              score(axts), symCount(axts), winSize, minScore, 
              as.character(resFiles),
              PACKAGE="CNEr")
  }
  CNE <- lapply(resFiles, 
               function(x){
                 res <- read.table(x, header=FALSE, sep="\t", as.is=TRUE)
               colnames(res) <- c("tName", "tStart", "tEnd", "qName", "qStart", 
                                 "qEnd", "strand", "score", "cigar")
               return(res)})
  names(CNE) <-  paste(minScore, winSize, sep="_")
  unlink(resFiles)
  return(CNE)
}


### -----------------------------------------------------------------
### Another main function for CNEs identification, 
### but it takes the axt files and bed files as input
### Not exported!
ceScanFile <- function(axtFiles, tFilterFile=NULL, qFilterFile=NULL, 
                       qSizes=NULL,
                       thresholds=c("49_50")){
  ## Here the returned tStart and qStart are 1-based coordinates. 
  ## Of course ends are also 1-based.
  if(!is.null(qFilterFile))
    if(is.null(qSizes) || !is(qSizes, "Seqinfo"))
      stop("qSizes must exist and be a Seqinfo object when qFilter exists")
  winSize <- as.integer(sapply(strsplit(thresholds, "_"), "[", 2))
  minScore <- as.integer(sapply(strsplit(thresholds, "_"), "[", 1))
  resFiles <- tempfile(pattern=paste(minScore, winSize, "ceScan", sep="-"), 
                      tmpdir=tempdir(), fileext="")
  .Call2("myCeScanFile", axtFiles, tFilterFile, qFilterFile, 
        as.character(seqnames(qSizes)), as.character(seqlengths(qSizes)),
        winSize, minScore,
        resFiles, PACKAGE="CNEr")
  CNE <- lapply(resFiles,
               function(x){
                 res <- read.table(x, header=FALSE, sep="\t", as.is=TRUE)
                 colnames(res) <- c("tName", "tStart", "tEnd", "qName", 
                                   "qStart", "qEnd", "strand", "score", "cigar")
                 return(res)})
  names(CNE) <-  paste(minScore, winSize, sep="_")
  unlink(resFiles)
  return(CNE)
}

### -----------------------------------------------------------------
### The S4 methods for ceScan
### Exported!
setMethod("ceScan", signature(axts="Axt", tFilter="GRanges", qFilter="GRanges",
                              qSizes="Seqinfo"),
          function(axts, tFilter, qFilter, qSizes, thresholds="49_50"){
            ceScanR(axts, tFilter=tFilter, qFilter=qFilter, 
                    qSizes=qSizes, thresholds=thresholds)
          }
          )
setMethod("ceScan", signature(axts="Axt", tFilter="missing", qFilter="GRanges",
                              qSizes="Seqinfo"),
          function(axts, tFilter, qFilter, qSizes, thresholds="49_50"){
            ceScanR(axts, tFilter=NULL, qFilter=qFilter,
                    qSizes=qSizes, thresholds=thresholds)
          }
          )
setMethod("ceScan", signature(axts="Axt", tFilter="missing", qFilter="missing",
                              qSizes="missing"),
          function(axts, tFilter, qFilter, qSizes, thresholds="49_50"){
            ceScanR(axts, tFilter=NULL, qFilter=NULL,
                    qSizes=NULL, thresholds=thresholds)
          }
          )
setMethod("ceScan", signature(axts="Axt", tFilter="GRanges", qFilter="missing",
                              qSizes="missing"),
          function(axts, tFilter, qFilter, qSizes, thresholds="49_50"){
            ceScanR(axts, tFilter=tFilter, qFilter=NULL,
                    qSizes=NULL, thresholds=thresholds)
          }
          )
setMethod("ceScan", signature(axts="character", tFilter="character", 
                              qFilter="character", qSizes="Seqinfo"),
          function(axts, tFilter, qFilter, qSizes, thresholds="49_50"){
            ceScanFile(axtFiles=axts, tFilterFile=tFilter, qFilterFile=qFilter,
                       qSizes=qSizes, thresholds=thresholds)
          }
          )
setMethod("ceScan", signature(axts="character", tFilter="missing",
                              qFilter="character", qSizes="Seqinfo"),
          function(axts, tFilter, qFilter, qSizes, thresholds="49_50"){
            ceScanFile(axtFiles=axts, tFilterFile=NULL, qFilterFile=qFilter,
                       qSizes=qSizes, thresholds=thresholds)
          }
          )
setMethod("ceScan", signature(axts="character", tFilter="missing",
                              qFilter="missing", qSizes="missing"),
          function(axts, tFilter, qFilter, qSizes, thresholds="49_50"){
            ceScanFile(axtFiles=axts, tFilterFile=NULL, qFilterFile=NULL,
                       qSizes=NULL, thresholds=thresholds)
          }
          )
setMethod("ceScan", signature(axts="character", tFilter="character",
                              qFilter="missing", qSizes="missing"),
          function(axts, tFilter, qFilter, qSizes, thresholds="49_50"){
            ceScanFile(axtFiles=axts, tFilterFile=tFilter, qFilterFile=NULL,
                       qSizes=NULL, thresholds=thresholds)
          }
          )

### -----------------------------------------------------------------
### Merge two side cnes
### Exported!
cneMerge <- function(cne1, cne2){
  # In this function, cne's start is 1-based coordinates. ends are 1-based too. 
  # Although in cne1 and cne2, query strand can be negative, 
  # but the coordinate is already on positive strand.
  ## first reverse the cne2's cigar
  ## cne2 = transform(cne2, cigar=chartr("DI", "ID", cigar))
  cne2$cigar <- chartr("DI", "ID", cne2$cigar)
  if(any(cne2$strand == "-")){
    #cne2[cne2$strand=="-", ] = transform(subset(cne2, strand=="-"), 
    #                                      cigar=reverseCigar(cigar))
    cne2[cne2$strand == "-", "cigar"] <- 
      reverseCigar(cne2[cne2$strand == "-", "cigar"])
  }
  colnames(cne2) <- c("qName", "qStart", "qEnd", "tName", 
                     "tStart", "tEnd", "strand", "score", "cigar")
  cne1T <- GRanges(seqnames=cne1$tName, 
                   ranges=IRanges(start=cne1$tStart, end=cne1$tEnd), 
                   strand="+")
  cne1Q <- GRanges(seqnames=cne1$qName, 
                   ranges=IRanges(start=cne1$qStart, end=cne1$qEnd), 
                   strand="+")
  cne2T <- GRanges(seqnames=cne2$tName, 
                   ranges=IRanges(start=cne2$tStart, end=cne2$tEnd), 
                   strand="+")
  cne2Q <- GRanges(seqnames=cne2$qName, 
                   ranges=IRanges(start=cne2$qStart, end=cne2$qEnd), 
                   strand="+")
  cneT <- c(cne1T, cne2T)
  cneQ <- c(cne1Q, cne2Q)
  # Here, I just removed the CNEs which are within another big CNEs. 
  # In very rare cases(1 in 100000), some cnes may just connect 
  # and need to merge them. Needs to be done in the future (perhaps not easy to be done in R).
  cneT_overlap <- findOverlaps(cneT, type="within", 
                              ignoreSelf=TRUE, ignoreRedundant=TRUE)
  #cneT_overlap1 = findOverlaps(cneT, type="equal", 
                              #ignoreSelf=TRUE, ignoreRedundant=TRUE)
  #cneT_overlap2 = findOverlaps(cneT, type="any", 
                              #ignoreSelf=TRUE, ignoreRedundant=TRUE)
  cneQ_overlap <- findOverlaps(cneQ, type="within", 
                              ignoreSelf=TRUE, ignoreRedundant=TRUE)
  #cneQ_overlap1 = findOverlaps(cneQ, type="equal", 
                              #ignoreSelf=TRUE, ignoreRedundant=TRUE)
  #cneQ_overlap2 = findOverlaps(cneQ, type="any", 
                              #ignoreSelf=TRUE, ignoreRedundant=TRUE)
  redundance <- IRanges::intersect(cneT_overlap, cneQ_overlap)
  #any_overlap = intersect(cneT_overlap2, cneQ_overlap2)
  #foo = setdiff(any_overlap, redundance)
  #paste(subjectHits(foo), queryHits(foo), sep=",") 
  # %in% paste(queryHits(redundance), subjectHits(redundance), sep=",")
  res <- rbind(cne1, cne2)[-queryHits(redundance), ] 
  # After the merge, we'd better name them as 1 and 2 
  # rather than the tName and qName. Use the names in mysql cne db.
  colnames(res) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", 
                    "strand", "similarity", "cigar")
  return(res)
}

#cutoffs1 = 4
#cutoffs2 = 8

blatCNE <- function(CNE, winSize, cutoffs1, cutoffs2, 
                    assembly1Twobit, assembly2Twobit, 
                    blatOptions=NULL, cutIdentity=90, 
                    tmpDir=tempdir(), blatBinary="blat"){
  # In this function, the input CNE's start and end are 1-based coordinates.
  blatOptionsALL <- list("DEF_BLAT_OPT_WSLO"=
                         "-tileSize=9 -minScore=24 -repMatch=16384",
                         "DEF_BLAT_OPT_WSMID"=
                         "-tileSize=10 -minScore=28 -repMatch=4096",
                         "DEF_BLAT_OPT_WSHI"=
                         "-tileSize=11 -minScore=30 -repMatch=1024")
  if(!is(winSize, "integer"))
    stop("winSize must be an integer!")
  if(is.null(blatOptions)){
    if(winSize > 45L)
      blatOptions <- blatOptionsALL[["DEF_BLAT_OPT_WSHI"]]
    else if(winSize > 35L)
      blatOptions <- blatOptionsALL[["DEF_BLAT_OPT_WSMID"]]
    else
      blatOptions <- blatOptionsALL[["DEF_BLAT_OPT_WSLO"]]
  }
  if(cutIdentity > 100  || cutIdentity < 0)
    stop("cutIdentity must be between 0 and 100!")
  if(!is(cutoffs1, "integer"))
    stop("cutoffs1 must be an integer!")
  if(!is(cutoffs2, "integer"))
    stop("cutoffs2 must be an integer!")
  
  .run_blat <- function(cne, cutIdentity, whichAssembly, 
                        assemblyTwobit, blatBinary, blatOptions, tmpDir){
    temp_cne <- tempfile(pattern="cne-", tmpdir=tmpDir)
    temp_psl <- tempfile(pattern="psl-", tmpdir=tmpDir)
    # For Blat, the start is 0-based and end is 1-based. 
    # So make cne's coordinates to comply with it.
    if(whichAssembly == 1){
      cne <- paste0(assemblyTwobit, ":", cne[,"chr1"], 
                    ":", cne[,"start1"]-1, "-", cne[,"end1"])
    }else{
      cne <- paste0(assemblyTwobit, ":", cne[,"chr2"], 
                    ":", cne[,"start2"]-1, "-", cne[,"end2"])
    }
    cne <- unique(cne)
    writeLines(cne, con=temp_cne)
    cmd <- paste0(blatBinary, " ", blatOptions, " ",
                  "-minIdentity=", cutIdentity,
                  " ", assemblyTwobit, " ", temp_cne, " ", temp_psl)
    my.system(cmd)
    unlink(temp_cne)
    return(temp_psl)
  }
  psl1Fn <- .run_blat(CNE, cutIdentity, 1, assembly1Twobit, 
                      blatBinary, blatOptions, tmpDir)
  psl2Fn <- .run_blat(CNE, cutIdentity, 2, assembly2Twobit, 
                      blatBinary, blatOptions, tmpDir)
  psl1 <- read.table(psl1Fn, header=FALSE, sep="\t", skip=5, as.is=TRUE)
  psl2 <- read.table(psl2Fn, header=FALSE, sep="\t", skip=5, as.is=TRUE)
  colnames(psl1) <- colnames(psl2) <- 
    c("matches", "misMatches", "repMatches", "nCount", 
      "qNumInsert", "qBaseInsert", "tNumInsert", "tBaseInsert", 
      "strand", "qName", "qSize", "qStart", "qEnd", "tName", 
      "tSize", "tStart", "tEnd", "blockCount", "blockSizes", 
      "qStarts", "tStarts")
  ##psl1 = subset(psl1, matches/qSize >= cutIdentity/100)
  ##for some reason, subset will cause error duing the Bioconductor build
  psl1 <- psl1[psl1$matches / psl1$qSize >= cutIdentity/100, ]
  ##psl2 = subset(psl2, matches/qSize >= cutIdentity/100)
  psl2 <- psl2[psl2$matches / psl2$qSize >= cutIdentity/100, ]
  psl1 <- table(psl1[ , "qName"])
  psl2 <- table(psl2[ , "qName"])
  CNEtNameIndex <- unname(psl1[paste0(CNE[ ,1], ":", CNE[ ,2]-1, "-", 
                                      CNE[ ,3])])
  CNEtNameIndex[is.na(CNEtNameIndex)] <- 0
  CNEqNameIndex <- unname(psl2[paste0(CNE[ ,4], ":", CNE[ ,5]-1, "-", 
                                      CNE[ ,6])])
  CNEqNameIndex[is.na(CNEqNameIndex)] <- 0
  CNE <- CNE[CNEtNameIndex <= cutoffs1 & CNEqNameIndex <= cutoffs2, ]
  return(CNE)
  # Here, the CNE's starts and ends are still 1-based.
}


ceScanOneStep <- function(axt1, filter1=NULL, sizes1, assembly1, twoBit1,
                          axt2, filter2=NULL, sizes2, assembly2, twoBit2,
                          thresholds=c("49_50"),
                          blatBinary="blat",
                          blatCutoff1, blatCutoff2
                          ){
  if(grepl("_", assembly1) || grepl("_", assembly2))
    stop("The assembly name must not contain \"_\"")
  if(assembly1 < assembly2){
    .ceScanSwap(axt1=axt1, filter1=filter1, sizes1=sizes1, assembly1=assembly1,
                twoBit1=twoBit1, axt2=axt2, filter2=filter2, sizes2=sizes2,
                assembly2=assembly2, twoBit2=twoBit2, thresholds=thresholds,
                blatBinary=blatBinary, blatCutoff1=blatCutoff1, 
                blatCutoff2=blatCutoff2)
  }else{
    .ceScanSwap(axt1=axt2, filter1=filter2, sizes1=sizes2, assembly1=assembly2,
                twoBit1=twoBit2, axt2=axt1, filter2=filter1,
                sizes2=sizes1, assembly2=assembly1, twoBit2=twoBit1,
                thresholds=thresholds, blatBinary=blatBinary,
                blatCutoff1=blatCutoff2, blatCutoff2=blatCutoff1)
  }
}

.ceScanSwap <- function(axt1, filter1=NULL, sizes1, assembly1, twoBit1,
                        axt2, filter2=NULL, sizes2, assembly2, twoBit2,
                        thresholds=c("49_50"),
                        blatBinary="blat",
                        blatCutoff1, blatCutoff2
                        ){
  ## In this function, we make sure "assembly1" is smaller than "assembly2".
  ## danRer7 is smaller than hg19.
  ## This is just for easier database storage table name: "danRer7_hg19_49_50"
  CNE1 <- ceScan(axt1, filter1, filter2, sizes2, thresholds)
  CNE2 <- ceScan(axt2, filter2, filter1, sizes1, thresholds)
  CNEMerged <- mapply(cneMerge, CNE1, CNE2, SIMPLIFY=FALSE)
  names(CNEMerged) = paste(assembly1, assembly2, thresholds, sep="_")
  CNEBlated <- list()
  for(i in 1:length(CNEMerged)){
    CNEBlated[[names(CNEMerged)[i]]] <- 
      blatCNE(CNEMerged[[i]], as.integer(sub(".+_.+_\\d+_", "", 
                                           names(CNEMerged)[i])),
              cutoffs1=blatCutoff1, cutoffs2=blatCutoff2,
              assembly1Twobit=twoBit1, assembly2Twobit=twoBit2,
              blatBinary=blatBinary)
  }
  ans <- CNE(assembly1=assembly1, assembly2=assembly2,
             thresholds=thresholds,
             CNE1=CNE1, CNE2=CNE2, CNEMerged=CNEMerged,
             CNERepeatsFiltered=CNEBlated,
             alignMethod=blatBinary)
  return(ans)
}



