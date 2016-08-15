### -----------------------------------------------------------------
### CNE length in power law
### Exported!
plotCNEWidth <- function(x, powerLawTest=FALSE, mc.cores=1L, ...){
  # x: GRangePairs
  if(!is(x, "GRangePairs")){
    stop("x must be a GRangePairs object!")
  }
  firstGRanges <- first(x)
  lastGRanges <- second(x)
  
  firstWidths <- width(firstGRanges)
  lastWidths <- width(lastGRanges)
  
  # Fitting a discrete power-law
  firstFit <- poweRlaw::displ$new(firstWidths)
  lastFit <- poweRlaw::displ$new(lastWidths)
  ## Estimate the x_min
  firstEst <- estimate_xmin(firstFit)
  lastEst <- estimate_xmin(lastFit)
  firstFit$setXmin(firstEst)
  lastFit$setXmin(lastEst)
  ## Estimate the scaling parameter alpha
  firstEst <- estimate_pars(firstFit)
  lastEst <- estimate_pars(lastFit)
  firstFit$setPars(firstEst)
  lastFit$setPars(lastEst)
  
  ## Prepare the xmin and pars test
  firstText <- paste("xmin:", firstFit$xmin, "\n",
                     "alpha:", format(firstFit$pars, digits = 3))
  lastText <- paste("xmin:", lastFit$xmin, "\n",
                     "alpha:", format(lastFit$pars, digits = 3))
  
  ## test power law distribution?
  if(powerLawTest){
    bs_p <- bootstrap_p(firstFit, no_of_sims=1000, threads=mc.cores)
    firstText <- paste(firstText, "\n", 
                       "pvalue:", format(bs_p$p, digits = 3))
    bs_p <- bootstrap_p(lastFit, no_of_sims=1000, threads=mc.cores)
    lastText <- paste(lastText, "\n", 
                       "pvalue:", format(bs_p$p, digits = 3))
  }
  
  ## Plot the distribution
  par(mfrow=c(1,2))
  poweRlaw::plot(firstFit, xlab="CNE width", ylab="CDF", main="first", pch=1, ...)
  poweRlaw::lines(firstFit, col=2, lwd=2)
  text(max(firstFit$dat), 0.7, firstText, adj = c( 1, 1 ))
  poweRlaw::plot(lastFit, xlab="CNE width", ylab="CDF", main="second", pch=1, ...)
  poweRlaw::lines(lastFit, col=2, lwd=2)
  text(max(lastFit$dat), 0.7, lastText, adj = c( 1, 1 ))
  
  invisible(list(first=firstFit, second=lastFit))
}

### -----------------------------------------------------------------
### CNE distribution:
### exported!
plotCNEDistribution <- function(x, chrs=NULL, chrScale=c("Mb", "Kb")){
  chrScale <- match.arg(chrScale)
  if(chrScale == "Mb"){
    chrScaleN <- 1e6
  }else{
    chrScaleN <- 1e3
  }
  ## x is a GRanges
  if(!is(x, "GRanges")){
    stop("`x` must be a `GRanges` object!")
  }
  
  if(seqlengthsNA(x)){
    stop("seqlengths must be provided in `x`!")
  }
  
  if(is.null(chrs)){
    ## Select the largest 6 chromosomes/scaffolds
    chrs <- names(sort(seqlengths(x), decreasing = TRUE))[1:6]
  }
  x <- x[seqnames(x) %in% chrs]
  x <- sort(x, ignore.strand=TRUE)
  dataToPlot <- data.frame(seqnames=seqnames(x),
                           x=round((start(x) + end(x))/2)/chrScaleN,
                           y=unlist(sapply(runLength(seqnames(x)), seq)))
  p <- ggplot(data=dataToPlot, aes_string(x="x",y="y")) +
    geom_point() + theme_bw() +
    facet_wrap(~seqnames, scales = "free") +
    xlab(paste0("Genomic location (", chrScale, ")")) + ylab("CNEs")
  p
}