

makeGRBs <- function(x, winSize=NULL, genes=NULL){
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

  # calculate the background percentage of coverage
  coveredAll <- sum(sum(cov)) / sum(as.numeric(seqlengths(x)))
  density <- runmean(cov, k=winSize,  endrule="constant")

  # slice the density into GRBs
  s <- slice(density, lower=coveredAll/2, includeLower=FALSE)
  clusterRanges <- GRanges(seqnames=rep(names(s), elementLengths(s)),
                           ranges=unlist(ranges(s)),
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

