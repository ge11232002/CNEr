### -----------------------------------------------------------------
### subAxt: get the subset of axt based on chr, start, end 
###         and from which sequence
### Exported!!
### -----------------------------------------------------------------
### subAxt method for character, integer, integer
setMethod("subAxt", signature(x="Axt", chr="character",
                              start="integer", end="integer"),
          function(x, chr, start, end, select=c("target", "query"),
                   qSize=NULL){
            searchGRanges <- GRanges(seqnames=chr,
                                     ranges=IRanges(start=start, end=end),
                                     strand="+")
            subAxt(x, searchGRanges, select=select, qSize=qSize)
          }
          )

### -----------------------------------------------------------------
### subAxt method for GRanges
setMethod("subAxt", signature(x="Axt", chr="GRanges",
                              start="missing", end="missing"),
          function(x, chr, start, end, select=c("target", "query"),
                   qSize=NULL){
            select <- match.arg(select)
            if(!is.null(qSize)){
              # If qSize is provided, use it to set the seqlengths of 
              # queryRanges
              seqlengths(second(x)) <- qSize[names(seqlengths(second(x)))]
            }
            # Fix the coordinates of queryRanges
            x <- fixCoordinates(x)
            
            strand(chr) <- "+"
            searchGRanges <- reduce(chr)
            ans <- .subAxtWhole(x, searchGRanges, select=select)
            
            # Restore the coordinates of queryRanges
            ans <- fixCoordinates(ans)
            return(ans)
          }
          )

### -----------------------------------------------------------------
### subAxt method for character, numeric, numeric
setMethod("subAxt", signature(x="Axt", chr="character",
                              start="numeric", end="numeric"),
          function(x, chr, start, end, select=c("target", "query"),
                   qSize=NULL){
            subAxt(x, chr, as.integer(start), as.integer(end), 
                   select=select, qSize=qSize)
          }
          )

### -----------------------------------------------------------------
### subAxt method for chr only
setMethod("subAxt", signature(x="Axt", chr="character",
                              start="missing", end="missing"),
          function(x, chr, start, end, select=c("target", "query"),
                   qSize=NULL){
            select <- match.arg(select)
            if(select == "target"){
              ans <- x[seqnames(targetRanges(x)) %in% chr]
            }else{
              ans <- x[seqnames(queryRanges(x)) %in% chr]
            }
            return(ans)
          }
          )

### -----------------------------------------------------------------
### .subAxtWhole: get the full axt alignment with searchGRanges
### Not exported!
.subAxtWhole <- function(x, searchGRanges, select=c("target", "query")){
  ## Here x is the Axt object with fixed coordinates.
  if(select == "target"){
    hitsAny <- findOverlaps(targetRanges(x),
                            searchGRanges, type="any",
                            select="all", ignore.strand=TRUE)
  }else if(select == "query"){
    hitsAny <- findOverlaps(queryRanges(x),
                            searchGRanges, type="any",
                            select="all", ignore.strand=TRUE)
  }
  indexAny <- queryHits(hitsAny)
  ans <- x[unique(indexAny)]
  return(ans)
}

### -----------------------------------------------------------------
### parallel subset of Axt alignment
### Exported!
psubAxt <- function(x, targetSearch, querySearch){
  stopifnot(is(x, "Axt"))
  if(!is(targetSearch, "GRanges") || !is(querySearch, "GRanges")){
    stop("targetSearch and querySearch must be `GRanges` object!")
  }
  if(length(targetSearch) != length(querySearch)){
    stop("targetSearch and querySearch must have same lengths.")
  }
  x <- fixCoordinates(x)
  
  hitsTarget <- findOverlaps(targetRanges(x), targetSearch, type="any",
                             select="all", ignore.strand=TRUE)
  hitsQuery <- findOverlaps(queryRanges(x), querySearch, type="any",
                            select="all", ignore.strand=TRUE)
  hits <- intersect(hitsTarget, hitsQuery)
  ans <- x[unique(queryHits(hits))]
  ans <- fixCoordinates(ans)
  return(ans)
}