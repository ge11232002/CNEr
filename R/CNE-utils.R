### -----------------------------------------------------------------
### CNE length power law
###
plotCNELength <- function(x){
  # x: GRangePairs
  firstGRanges <- first(x)
  lastGRanges <- last(x)
  
  firstWidths <- width(firstGRanges)
  firstWidthsLog <- log10(firstWidths)
  x <- seq(1.5, 3.1, by=0.05)
  y <- log10(sapply(x, function(x){sum(firstWidthsLog >= x)}))
 # plot(x,y)
}