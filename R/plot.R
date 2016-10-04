### -----------------------------------------------------------------
### Fetch the CNE coordinates from SQL and compute the densities
### Exported!
CNEDensity <- function(dbName, tableName, chr, start, end, 
                       whichAssembly=c("first", "second"),
                       windowSize=300, minLength=NULL){
  .CNEDensityInternal(dbName=dbName, tableName=tableName,
                      whichAssembly=whichAssembly,
                      chr=chr, start=start, end=end,
                      windowSize=windowSize, minLength=minLength)
}

.CNEDensityInternal <- function(dbName, tableName,
                                whichAssembly=c("first","second"),
                                chr, start, end, windowSize, 
                                minLength=NULL){
  nrGraphs <- 1
  CNEstart <- start
  CNEend <- end
  stopifnot(length(CNEstart) == 1L)
  stopifnot(length(CNEend) == 1L)
  # This is the pipeline of doing the density plot
  # The windowSize is in kb.
  whichAssembly <- match.arg(whichAssembly)
  windowSize <- as.integer(windowSize) * 1000
  # make things easier
  if(windowSize %% 2 == 0)
    windowSize <- windowSize - 1L
  context_start <- as.integer(max(CNEstart - (windowSize-1L)/2, 1))
  context_end <- as.integer(CNEend + (windowSize-1)/2)
  rangesPair <- readCNERangesFromSQLite(dbName, tableName, chr,
                                    context_start, context_end, 
                                    whichAssembly, minLength)
  ## When no CNEs are returned
  if(length(rangesPair) == 0L){
    ans <- GRanges(seqnames=chr,
                   ranges=IRanges(start=context_start,
                                  end=context_end),
                   strand="*",
                   score=0)
    return(ans)
  }
  
  if(whichAssembly == "first"){
    ranges <- first(rangesPair)
  }else if(whichAssembly == "second"){
    ranges <- second(rangesPair)
  }
  # Implement get_cne_ranges_in_region_partitioned_by_other_chr later!!!
  ranges <- reduce(ranges)
  covAll <- coverage(ranges, width=context_end)
  runMeanAll <- runmean(covAll, k=windowSize, "constant")
  ans <- as(runMeanAll, "GRanges")
  ans$score <- ans$score * 100
  return(ans)
}