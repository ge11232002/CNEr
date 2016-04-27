### -----------------------------------------------------------------
### make GRBs from CNEs
### Exported!
makeGRBs <- function(x, winSize=NULL, genes=NULL, ratio=0.5){
  if(!is(x, "GRangesList")){
    stop("The input `x` must be a `GRangesList` object!")
  }
  
  # winSize in kb
  x <- unlist(x)
  x <- reduce(x)
  cov <- coverage(x)

  # Guess winSize from genome size if NULL
  if(is.null(winSize)){
    ## Based on 300kb for human
    winSize <- seqlengths(x) / 3e6 * 300
  }
  winSize <- winSize * 1e3
  
  if(!is.null(genes) && !is(genes, "GRanges")){
    stop("The `genes` must be a `GRanges` object!")
  }

  # calculate the background percentage of coverage
  coveredAll <- sum(sum(cov)) / sum(as.numeric(seqlengths(x)))
  density <- runmean(cov, k=winSize,  endrule="constant")

  # slice the density into GRBs
  slicedDensities <- slice(density, lower=coveredAll*ratio, includeLower=FALSE)
  clusterRanges <- GRanges(seqnames=rep(names(slicedDensities), 
                                        elementNROWS(slicedDensities)),
                           ranges=unlist(ranges(slicedDensities)),
                           strand="+",
                           seqinfo=seqinfo(x))

  # remove GRBs (which do not encompass any gene)
  if(!is.null(genes)){
    hits <- findOverlaps(genes, clusterRanges, type="within", select="all",
                         ignore.strand=TRUE)
    indexKeep <- unique(subjectHits(hits))
    clusterRanges <- clusterRanges[indexKeep]
  }
  return(clusterRanges)
}

