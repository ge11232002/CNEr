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
            if(select == "query"){
              if(is.null(qSize)){
                qSize <- seqlengths(queryRanges(x))
              }
              if(!is(qSize, "integer")){
                stop("qSize must be an integer object.")
              }
              if(!all(as.character(seqnames(chr)) %in% names(qSize))){
                stop("All the chromosomes in `x` must exist in `qSize`.")
              }
            }
            strand(chr) <- "+"
            searchGRanges <- reduce(chr)
            .subAxtWhole(x, searchGRanges, select=select, qSize=qSize)
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
              ans <- x[as.character(seqnames(targetRanges(x))) %in% chr]
            }else{
              ans <- x[as.character(seqnames(queryRanges(x))) %in% chr]
            }
            return(ans)
          }
          )

### -----------------------------------------------------------------
### .subAxtWhole: get the full axt alignment with searchGRanges
### Not exported!
.subAxtWhole <- function(x, searchGRanges, select=c("target", "query"),
                         qSize=NULL){
  if(select == "target"){
    # It's easy because the strand on target reference is always positive.
    hitsAny <- findOverlaps(targetRanges(x),
                            searchGRanges, type="any",
                            select="all", ignore.strand=TRUE)
    indexAny <- queryHits(hitsAny)
    ans <- x[unique(indexAny)]
  }else if(select == "query"){
    start <- start(searchGRanges)
    end <- end(searchGRanges)
    # first search Axts on positive strand
    ## searchGRanges has posive strands inside.
    hitsPositiveAny <- findOverlaps(queryRanges(x),
                                    searchGRanges, type="any",
                                    select="all", ignore.strand=FALSE)
    indexPositiveAny <- queryHits(hitsPositiveAny)
    
    qSize <- qSize[as.character(seqnames(searchGRanges))]
    
    # then search Axts on negative strand.
    ## we need to prepare the searchGRanges on negative strand.
    searchGRangesNegative <- GRanges(seqnames=seqnames(searchGRanges),
                                     ranges=IRanges(start=qSize-end+1,
                                                    end=qSize-start+1),
                                     strand="-")
    searchGRangesNegative <- reduce(searchGRangesNegative)
    hitsNegativeAny <- findOverlaps(queryRanges(x),
                                    searchGRangesNegative, type="any",
                                    select="all", ignore.strand=FALSE)
    indexNegativeAny <- queryHits(hitsNegativeAny)
    ans <- x[unique(c(indexPositiveAny, indexNegativeAny))]
  }
  return(ans)
}