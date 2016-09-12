### -----------------------------------------------------------------
### The main function for scanning the axts and get the CNEs
### Not exported!
ceScanR <- function(axts, tFilter=NULL, qFilter=NULL, tSizes, qSizes, 
                    thresholds=c("49_50")){
  winSize <- as.integer(sapply(strsplit(thresholds, "_"), "[", 2))
  minScore <- as.integer(sapply(strsplit(thresholds, "_"), "[", 1))
  resFiles <- tempfile(pattern=paste(minScore, winSize, "ceScan", sep="-"), 
                       tmpdir=tempdir(), fileext="")
  ## How stupid I am...
  if(is.null(tFilter) && is.null(qFilter)){
    .Call2("myCeScan", NULL, NULL, NULL,
           NULL, NULL, NULL,
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
  }else{
    .Call2("myCeScan", 
           as.character(seqnames(tFilter)), start(tFilter), end(tFilter),
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
                  res <- try(read.table(x, header=FALSE, sep="\t", as.is=TRUE),
                             silent=TRUE)
                  if(class(res) == "try-error"){
                    return(GRangePairs(first=GRanges(seqinfo=tSizes), 
                                       second=GRanges(seqinfo=qSizes))
                           )
                  }
                  colnames(res) <- c("tName", "tStart", "tEnd", 
                                     "qName", "qStart", "qEnd",
                                     "strand", "score", "cigar")
                  ans <- GRangePairs(first=GRanges(seqnames=res$tName,
                                             ranges=IRanges(start=res$tStart,
                                                            end=res$tEnd),
                                             strand="+",
                                             seqinfo=tSizes),
                                     second=GRanges(seqnames=res$qName,
                                              ranges=IRanges(start=res$qStart,
                                                             end=res$qEnd),
                                              strand=res$strand,
                                              seqinfo=qSizes),
                                     score=res$score,
                                     cigar=res$cigar
                                     )
                  return(ans)
                  }
                )
  names(CNE) <-  paste(minScore, winSize, sep="_")
  unlink(resFiles)
  return(CNE)
}

### -----------------------------------------------------------------
### The function for ceScan: one way for axt object, two way for CNE
### Exported!
setGeneric("ceScan", function(x, tFilter=NULL, qFilter=NULL,
                              tSizes=NULL, qSizes=NULL,
                              window=50L, identity=50L)
  standardGeneric("ceScan"))

setMethod("ceScan", "Axt", function(x, tFilter=NULL, qFilter=NULL,
                                    tSizes=NULL, qSizes=NULL,
                                    window=50L, identity=50L){
  ceScanAxt(x, tFilter=tFilter, qFilter=qFilter, tSizes=tSizes, qSizes=qSizes, 
            window=window, identity=identity)
})

setMethod("ceScan", "CNE", function(x, tFilter=NULL, qFilter=NULL,
                                    tSizes=NULL, qSizes=NULL,
                                    window=50L, identity=50L){
  axt12 <- readAxt(x@axt12Fn)
  axt21 <- readAxt(x@axt21Fn)
  cne12 <- ceScan(axt12, tFilter, qFilter,
                  tSizes=seqinfo(TwoBitFile(x@assembly1Fn)),
                  qSizes=seqinfo(TwoBitFile(x@assembly2Fn)),
                  window=window, identity=identity)
  cne21 <- ceScan(axt21, qFilter, tFilter,
                  tSizes=seqinfo(TwoBitFile(x@assembly2Fn)),
                  qSizes=seqinfo(TwoBitFile(x@assembly1Fn)),
                  window=window, identity=identity)
  ans <- list()
  for(i in 1:length(cne12)){
    ans[[names(cne12)[i]]] <- BiocGenerics:::replaceSlots(x, 
                                                          CNE12=cne12[[i]],
                                                          CNE21=cne21[[i]],
                                                          window=window[i],
                                                          identity=identity[i])
  }
  return(ans)
})

ceScanAxt <- function(axts, tFilter=NULL, qFilter=NULL,
                      tSizes=NULL, qSizes=NULL,
                      window=50L, identity=50L){
  # Prepare tSizes
  if(is.null(tSizes)){
    ## tSizes is NULL
    tSizes <- seqinfo(targetRanges(axts))
    if(any(is.na(seqlengths(tSizes)))){
      stop("tSizes must be provided or seqinfo is available in Axt object!")
    }
  }else if(!is(tSizes, "Seqinfo")){
    ## tSizes is integer vector
    tSizes <- Seqinfo(seqnames=names(tSizes), seqlengths=tSizes)
  }
  ## Check tFilter, tFilter's seqnames must %in% tSizes
  if(!is.null(tFilter)){
    if(!is(tFilter, "GRanges")){
      stop("tFilter must be NULL or a GRanges object!")
    }else{
      if(!all(seqnames(tFilter) %in% seqnames(tSizes))){
        stop("All the chromosomes in tFilter must exist in tSizes!")
      }
    }
  }
  ## Check axts target chromosomes
  if(!all(seqnames(targetRanges(axts)) %in% seqnames(tSizes))){
    stop("All the chromosomes in targetRanges of Axt object must exist in tSizes!")
  }
  
  # Prepare qSizes
  if(is.null(qSizes)){
    ## qSizes is NULL
    qSizes <- seqinfo(queryRanges(axts))
    if(any(is.na(seqlengths(qSizes)))){
      stop("qSizes must be provided or seqinfo is available in Axt object!")
    }
  }else if(!is(qSizes, "Seqinfo")){
    ## qSizes is integer vector
    qSizes <- Seqinfo(seqnames=names(qSizes), seqlengths=qSizes)
  }
  ## Check qFilter, qFilter's seqnames must %in% qSizes
  if(!is.null(qFilter)){
    if(!is(qFilter, "GRanges")){
      stop("qFilter must be NULL or a GRanges object!")
    }else{
      if(!all(seqnames(qFilter) %in% seqnames(qSizes))){
        stop("All the chromosomes in qFilter must exist in qSizes!")
      }
    }
  }
  ## Check axts query chromosomes
  if(!all(seqnames(queryRanges(axts)) %in% seqnames(qSizes))){
    stop("All the chromosomes in queryRanges of Axt object must exist in qSizes!")
  }
  
  # Check identity and window size
  if(any(window < identity)){
    stop("The scanning window size must be equal or larger than identity!")
  }
  if(length(identity) != length(window)){
    stop("The length of identity and window are different!")
  }
  thresholds <- paste(identity, window, sep="_")
  ans <- ceScanR(axts, tFilter=tFilter, qFilter=qFilter, 
                 tSizes=tSizes, qSizes=qSizes, thresholds=thresholds)
}

### -----------------------------------------------------------------
### Merge two side cnes (GRangePairs object), mcols will be discarded.
### strand information is no longer important.
### Exported!
setGeneric("cneMerge", function(cne12, cne21) standardGeneric("cneMerge"))

setMethod("cneMerge", signature(cne12="CNE", cne21="missing"),
          function(cne12, cne21){
            cneMerged <- cneMergeGRangePairs(cne12@CNE12, cne12@CNE21)
            BiocGenerics:::replaceSlots(cne12, CNEMerged=cneMerged)
          })

setMethod("cneMerge", signature(cne12="GRangePairs", cne21="GRangePairs"),
          function(cne12, cne21){
            cneMergeGRangePairs(cne12, cne21)
          })

cneMergeGRangePairs <- function(cne12, cne21){
  if(!is(cne12, "GRangePairs") || !is(cne21, "GRangePairs")){
    stop("cne12 and cne21 must be a GRangePairs object!")
  }
  strand(cne12@first) <- strand(cne12@second) <- strand(cne21@first) <- 
    strand(cne21@second) <- "+"
  cne <- c(cne12, swap(cne21))
  
  # 1. using reduce: the problem: it won't deal with {1,2,3}, {1,2} case
  # firstReduce <- reduce(first(cne), with.revmap=TRUE)
  # lastReduce <- reduce(last(cne), with.revmap=TRUE)
  # overlapFirstLast <- IntegerList(intersect(firstReduce$revmap,
  #                                           lastReduce$revmap))
  # 
  # ## First deal with the merged ranges
  # overlapFirstLast <- overlapFirstLast[elementNROWS(overlapFirstLast) > 1L]
  # firstIndex <- match(paste(overlapFirstLast, collapse="-"),
  #                     paste(firstReduce$revmap, collapse="-"))
  # lastIndex <- match(paste(overlapFirstLast, collapse="-"),
  #                    paste(lastReduce$revmap, collapse="-"))
  # ansFirst <- firstReduce[firstIndex]
  # mcols(ansFirst) <- NULL
  # ansLast <- lastReduce[lastIndex]
  # mcols(ansLast) <- NULL
  # 
  # ## Then deal with the unmerged ranges
  # ansFirst <- c(ansFirst, first(cne)[-sort(unlist(overlapFirstLast))])
  # ansLast <- c(ansLast, last(cne)[-sort(unlist(overlapFirstLast))])
  # return(GRangePairs(first=ansFirst, last=ansLast))
  
  # 2. using the findOverlaps:
  firstHits <- findOverlaps(first(cne), type="within",
                            drop.self=TRUE, drop.redundant=TRUE)
  lastHist <- findOverlaps(second(cne), type="within",
                           drop.self=TRUE, drop.redundant=TRUE)
  redundance <- IRanges::intersect(firstHits, lastHist)
  if(length(redundance) == 0L){
    return(cne)
  }else{
    return(cne[-queryHits(redundance)])
  }
}

### -----------------------------------------------------------------
### blatCNE: 
### Exported!
blatCNE <- function(cne, blatOptions=NULL, cutIdentity=90){
  if(!is(cne, "CNE")){
    stop("CNE must be a CNE class object!")
  }
  if(!file.exists(cne@assembly1Fn)){
    stop("The assembly1 twoBit file must exit!")
  }
  if(!file.exists(cne@assembly2Fn)){
    stop("The assembly2 twoBit file must exit!")
  }
  blatOptionsALL <- list("DEF_BLAT_OPT_WSLO"=
                         "-tileSize=9 -minScore=24 -repMatch=16384",
                         "DEF_BLAT_OPT_WSMID"=
                         "-tileSize=10 -minScore=28 -repMatch=4096",
                         "DEF_BLAT_OPT_WSHI"=
                         "-tileSize=11 -minScore=30 -repMatch=1024")
  winSize <- cne@window
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
  cutoffs1 <- cne@cutoffs1
  cutoffs2 <- cne@cutoffs2
  
  .run_blat <- function(cne, cutIdentity, whichAssembly=c("first", "last"),
                        blatOptions){
    whichAssembly <- match.arg(whichAssembly)
    temp_cne <- tempfile(pattern="cne-")
    temp_psl <- tempfile(pattern="psl-")
    # For Blat, the start is 0-based and end is 1-based. 
    # So make cne's coordinates to comply with it.
    if(whichAssembly == "first"){
      assemblyTwobit <- cne@assembly1Fn
      cneDataFrame <- paste0(assemblyTwobit, ":", 
                    as.character(seqnames(first(cne@CNEMerged))),
                    ":", format(start(first(cne@CNEMerged))-1,
                                trim=TRUE, scientific=FALSE),
                    "-", format(end(first(cne@CNEMerged)),
                                trim=TRUE, scientific=FALSE))
    }else{
      assemblyTwobit <- cne@assembly2Fn
      cneDataFrame <- paste0(assemblyTwobit, ":", 
                    as.character(seqnames(second(cne@CNEMerged))), 
                    ":", format(start(second(cne@CNEMerged))-1,
                                trim=TRUE, scientific=FALSE),
                    "-", format(end(second(cne@CNEMerged)),
                                trim=TRUE, scientific=FALSE))
    }
    cneDataFrame <- unique(cneDataFrame)
    writeLines(cneDataFrame, con=temp_cne)
    cmd <- paste0(cne@aligner, " ", blatOptions, " ",
                  "-minIdentity=", cutIdentity,
                  " ", assemblyTwobit, " ", temp_cne, " ", temp_psl)
    my.system(cmd)
    unlink(temp_cne)
    return(temp_psl)
  }
  psl1Fn <- .run_blat(cne, cutIdentity, "first", blatOptions)
  psl2Fn <- .run_blat(cne, cutIdentity, "last", blatOptions)
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
  CNEtNameIndex <- unname(psl1[paste0(as.character(
                                        seqnames(first(cne@CNEMerged))),
                                      ":", format(start(first(cne@CNEMerged))-1,
                                                  trim=TRUE, scientific=FALSE),
                                      "-", format(end(first(cne@CNEMerged)),
                                                  trim=TRUE, scientific=FALSE))
                               ]
                          )
  CNEtNameIndex[is.na(CNEtNameIndex)] <- 0
  CNEqNameIndex <- unname(psl2[paste0(as.character(
                                        seqnames(second(cne@CNEMerged))), 
                                      ":", format(start(second(cne@CNEMerged))-1,
                                                  trim=TRUE, scientific=FALSE),
                                      "-", format(end(second(cne@CNEMerged)),
                                                  trim=TRUE, scientific=FALSE))
                               ]
                          )
  CNEqNameIndex[is.na(CNEqNameIndex)] <- 0
  cneFinal <- cne@CNEMerged[as.integer(CNEtNameIndex) <= cutoffs1 &
                            as.integer(CNEqNameIndex) <= cutoffs2]
  BiocGenerics:::replaceSlots(cne, CNEFinal=cneFinal)
}

# ceScanOneStep <- function(axt1, filter1=NULL, sizes1, assembly1, twoBit1,
#                           axt2, filter2=NULL, sizes2, assembly2, twoBit2,
#                           thresholds=c("49_50"),
#                           blatBinary="blat",
#                           blatCutoff1, blatCutoff2
#                           ){
#   if(grepl("_", assembly1) || grepl("_", assembly2))
#     stop("The assembly name must not contain \"_\"")
#   if(assembly1 < assembly2){
#     .ceScanSwap(axt1=axt1, filter1=filter1, sizes1=sizes1, assembly1=assembly1,
#                 twoBit1=twoBit1, axt2=axt2, filter2=filter2, sizes2=sizes2,
#                 assembly2=assembly2, twoBit2=twoBit2, thresholds=thresholds,
#                 blatBinary=blatBinary, blatCutoff1=blatCutoff1, 
#                 blatCutoff2=blatCutoff2)
#   }else{
#     .ceScanSwap(axt1=axt2, filter1=filter2, sizes1=sizes2, assembly1=assembly2,
#                 twoBit1=twoBit2, axt2=axt1, filter2=filter1,
#                 sizes2=sizes1, assembly2=assembly1, twoBit2=twoBit1,
#                 thresholds=thresholds, blatBinary=blatBinary,
#                 blatCutoff1=blatCutoff2, blatCutoff2=blatCutoff1)
#   }
# }
# 
# .ceScanSwap <- function(axt1, filter1=NULL, sizes1, assembly1, twoBit1,
#                         axt2, filter2=NULL, sizes2, assembly2, twoBit2,
#                         thresholds=c("49_50"),
#                         blatBinary="blat",
#                         blatCutoff1, blatCutoff2
#                         ){
#   ## In this function, we make sure "assembly1" is smaller than "assembly2".
#   ## danRer7 is smaller than hg19.
#   ## This is just for easier database storage table name: "danRer7_hg19_49_50"
#   CNE1 <- ceScan(axt1, filter1, filter2, sizes2, thresholds)
#   CNE2 <- ceScan(axt2, filter2, filter1, sizes1, thresholds)
#   CNEMerged <- mapply(cneMerge, CNE1, CNE2, SIMPLIFY=FALSE)
#   names(CNEMerged) = paste(assembly1, assembly2, thresholds, sep="_")
#   CNEBlated <- list()
#   for(i in 1:length(CNEMerged)){
#     CNEBlated[[names(CNEMerged)[i]]] <- 
#       blatCNE(CNEMerged[[i]], as.integer(sub(".+_.+_\\d+_", "", 
#                                            names(CNEMerged)[i])),
#               cutoffs1=blatCutoff1, cutoffs2=blatCutoff2,
#               assembly1Twobit=twoBit1, assembly2Twobit=twoBit2,
#               blatBinary=blatBinary)
#   }
#   ans <- CNE(assembly1=assembly1, assembly2=assembly2,
#              thresholds=thresholds,
#              CNE1=CNE1, CNE2=CNE2, CNEMerged=CNEMerged,
#              CNERepeatsFiltered=CNEBlated,
#              alignMethod=blatBinary)
#   return(ans)
# }
