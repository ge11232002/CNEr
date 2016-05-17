### -----------------------------------------------------------------
### make GRBs from CNEs
### Exported!
makeGRBs <- function(x, winSize=NULL, genes=NULL, ratio=0.5){
  if(!is(x, "GRangesList")){
    stop("The input `x` must be a `GRangesList` object!")
  }
  xList <- x
  if(is.null(names(xList))){
    names(xList) <- as.character(1:length(xList))
  }
  x <- unlist(x)
  x <- reduce(x)
  cov <- coverage(x)

  # Guess winSize from genome size if NULL
  # winSize in kb
  if(is.null(winSize)){
    ## Based on 300kb for human
    winSize <- seqlengths(x) / 3e6 * 300
  }
  winSize <- winSize * 1e3
  
  if(!is.null(genes) && !is(genes, "GRanges")){
    stop("The `genes` must be a `GRanges` object!")
  }

  # calculate the background percentage of coverage
  totalGenomeSize <- 
    sum(as.numeric(seqlengths(x)[as.character(unique(seqnames(x)))]))
  if(is.na(totalGenomeSize)){
    stop("seqlengths must be provided in input x!")
  }
  coveredAll <- sum(sum(cov)) / totalGenomeSize
  density <- runmean(cov, k=winSize,  endrule="constant")

  # slice the density into GRBs
  slicedDensities <- slice(density, lower=coveredAll*ratio, includeLower=FALSE)
  clusterRanges <- GRanges(seqnames=rep(names(slicedDensities), 
                                        elementNROWS(slicedDensities)),
                           ranges=unlist(ranges(slicedDensities)),
                           strand="+",
                           seqinfo=seqinfo(x))

  # shrink the GRBs with actual CNE locations
  hits <- findOverlaps(x, clusterRanges, type="within", select="all",
                       ignore.strand=TRUE)
  mergedGRBCNEsList <- split(x[queryHits(hits)], subjectHits(hits))
  starts <- min(start(mergedGRBCNEsList))
  ends <- max(end(mergedGRBCNEsList))
  seqnames <- sapply(seqnames(mergedGRBCNEsList), runValue)
  GRBMergedClean <- GRanges(seqnames=seqnames,
                            ranges=IRanges(start=starts,
                                           end=ends),
                            strand="+",
                            seqinfo=seqinfo(x))
  clusterRanges <- GRBMergedClean
  
  # remove GRBs (which do not encompass any gene)
  if(!is.null(genes)){
    hits <- findOverlaps(genes, clusterRanges, type="within", select="all",
                         ignore.strand=TRUE)
    indexKeep <- unique(subjectHits(hits))
    clusterRanges <- clusterRanges[indexKeep]
  }
  
  # Count the number of CNEs within the GRBs
  for(i in 1:length(xList)){
    hitsCNEs <- findOverlaps(xList[[i]], clusterRanges,
                             ignore.strand=TRUE, type="within")
    cnes <- sapply(split(queryHits(hitsCNEs), subjectHits(hitsCNEs)), length)
    mcols(clusterRanges)[[names(xList)[i]]] <- 0L
    mcols(clusterRanges)[[names(xList)[i]]][as.integer(names(cnes))] <- cnes
  }
  
  return(clusterRanges)
}

