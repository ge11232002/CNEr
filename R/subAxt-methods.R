### -----------------------------------------------------------------
### get the subset of axt based on chr, start, end and from which sequence
### Exported!!

setMethod("subAxt", signature(x="Axt", chr="character",
                              start="integer", end="integer"),
          function(x, chr, start, end, select=c("target", "query"), qSize=NULL){
            .subAxtMultiple(x, chr=chr, start=start, end=end,
                            select=select, qSize=qSize)
          }
          )

setMethod("subAxt", signature(x="Axt", chr="character",
                              start="numeric", end="numeric"),
          function(x, chr, start, end, select=c("target", "query"), qSize=NULL){
            subAxt(x, chr, as.integer(start), as.integer(end), 
                   select=select, qSize=qSize)
          }
          )

setMethod("subAxt", signature(x="Axt", chr="character",
                              start="missing", end="missing"),
          function(x, chr, start, end, select=c("target", "query"), qSize=NULL){
            select <- match.arg(select)
            if(select == "target"){
              ans <- x[seqnames(targetRanges(x)) %in% chr]
            }else{
              ans <- x[seqnames(queryRanges(x)) %in% chr]
            }
            return(ans)
          }
          )

## This is to fetch the Axts within the specific chrs, starts, ends
## based on target sequences.
.subAxtFull <- function(x, chr, start, end,
                        select=c("target", "query"),
                        #type=c("any", "within"),
                        qSize=NULL){
  #type <- match.arg(type)
  select <- match.arg(select)
  if(select == "query"){
    if(is.null(qSize))
      stop("When selecting on query alignments,
           qSize must be provided.")
           if(!is(qSize, "integer"))
             stop("qSize must be an integer object.")
  }
  if(length(chr) != 1L || length(start) != 1L || length(end) != 1L){
    stop("Currently we only support chr, start, end in length of 1.")
  }
  start = as.integer(start)
  end = as.integer(end)

  if(select == "target"){
    searchGRanges <- GRanges(seqnames=chr,
                             ranges=IRanges(start=start, end=end),
                             strand="+")
    ## First we find the axts totally within the coordinates
    indexWithin <- queryHits(findOverlaps(targetRanges(x),
                                          searchGRanges, type="within",
                                          select="all"))
    indexAny <- queryHits(findOverlaps(targetRanges(x),
                                       searchGRanges, type="any",
                                       select="all"))
    ## Specially take care of the axt just partially overlapped with the range
    indexPartial <- setdiff(indexAny, indexWithin)
    newAxts = list()
    if(length(indexPartial) > 0L){
      newStarts <- pmax(start(targetRanges(x)[indexPartial]), start)
      newEnds <- pmin(end(targetRanges(x)[indexPartial]), end)
      for(i in 1:length(indexPartial)){
        newAxts[[i]] <- subAln(x[indexPartial[i]], newStarts[i], newEnds[i],
                               select="target")
      }
    }
    ans <- c(x[indexWithin], do.call(c, newAxts))
    return(ans)
  }else if(select == "query"){
    ## first search Axts on positive strand
    searchGRanges <- GRanges(seqnames=chr,
                             ranges=IRanges(start=start, end=end),
                             strand="+")
    indexPositiveWithin <- queryHits(findOverlaps(queryRanges(x),
                                                  searchGRanges, type="within",
                                                  select="all"))
    indexPositiveAny <- queryHits(findOverlaps(queryRanges(x),
                                               searchGRanges, type="any",
                                               select="all"))
    indexPartial <- setdiff(indexPositiveAny, indexPositiveWithin)
    newAxts = list()
    if(length(indexPartial) > 0L){
      newStarts <- pmax(start(queryRanges(x))[indexPartial], start)
      newEnds <- pmin(end(queryRanges(x))[indexPartial], end)
      for(i in 1:length(indexPartial)){
        newAxts[[i]] <- subAln(x[indexPartial[i]], newStarts[i], newEnds[i],
                               select="query")
      }
    }
    # then search Axts on negative strand.
    # we need to prepare the searchGRanges on negative strand.
    searchGRanges <- GRanges(seqnames=chr,
                             ranges=IRanges(start=qSize-end+1,
                                            end=qSize-start+1),
                             strand="-")
    indexNegativeWithin <- queryHits(findOverlaps(queryRanges(x),
                                                  searchGRanges, type="within",
                                                  select="all"))
    indexNegativeAny <- queryHits(findOverlaps(queryRanges(x),
                                               searchGRanges, type="any",
                                               select="all"))
    indexPartial <- setdiff(indexNegativeAny, indexNegativeWithin)
    newAxts2 = list()
    if(length(indexPartial) > 0L){
      newStarts <- pmax(start(queryRanges(x))[indexPartial], qSize-end+1)
      newEnds <- pmin(end(queryRanges(x))[indexPartial], qSize-start+1)
      for(i in 1:length(indexPartial)){
        newAxts2[[i]] <- subAln(x[indexPartial[i]], newStarts[i], newEnds[i],
                                select="query")
      }
    }
    index <- sort(c(indexPositiveWithin, indexNegativeWithin))
    ans <- c(x[index], do.call(c, newAxts), do.call(c, newAxts2))
    return(ans)
  }
}

.subAxtMultiple <- function(x, chr, start, end,
                            select=c("target", "query"),
                            qSize=NULL){
  ## We will allow multiple chr, start, end
  select <- match.arg(select)
  if(select == "query"){
    if(is.null(qSize)){
      stop("When selecting on query alignments,
           qSize must be provided.")
    }
    if(!is(qSize, "integer")){
      stop("qSize must be an integer object.")
    }
  }
  searchGRanges <- GRanges(seqnames=chr,
                           ranges=IRanges(start=start, end=end),
                           strand="+")
  searchGRanges <- reduce(searchGRanges)
  start <- start(searchGRanges)
  end <- end(searchGRanges)
  if(select == "target"){
    ## First we find the axts totally within the coordinates
    hitsWithin <- findOverlaps(targetRanges(x),
                               searchGRanges, type="within",
                               select="all")
    indexWithin <- queryHits(hitsWithin)
    hitsAny <- findOverlaps(targetRanges(x),
                            searchGRanges, type="any",
                            select="all")
    indexAny <- queryHits(hitsAny)
    ## Specially take care of the axt just partially overlapped with the range
    hitsPartial <- setdiff(hitsAny, hitsWithin)
    newAxts = list()
    if(length(hitsPartial) > 0L){
      newStarts <- pmax(start(targetRanges(x)[queryHits(hitsPartial)]), 
                        start(searchGRanges)[subjectHits(hitsPartial)])
      newEnds <- pmin(end(targetRanges(x)[queryHits(hitsPartial)]), 
                      end(searchGRanges)[subjectHits(hitsPartial)])
      for(i in 1:length(hitsPartial)){
        newAxts[[i]] <- subAln(x[queryHits(hitsPartial[i])], 
                               newStarts[i], newEnds[i],
                               select="target")
      }
    }
    ans <- c(x[indexWithin], do.call(c, newAxts))
    return(ans)
  }else{
    ## first search Axts on positive strand
    hitsPositiveWithin <- findOverlaps(queryRanges(x),
                                       searchGRanges, type="within",
                                       select="all")
    indexPositiveWithin <- queryHits(hitsPositiveWithin)
    hitsPositiveAny <- findOverlaps(queryRanges(x),
                                    searchGRanges, type="any",
                                    select="all")
    indexPositiveAny <- queryHits(hitsPositiveAny)
    hitsPartial <- setdiff(hitsPositiveAny, hitsPositiveWithin)
    newAxts = list()
    if(length(hitsPartial) > 0L){
      newStarts <- pmax(start(queryRanges(x))[queryHits(hitsPartial)], 
                        start(searchGRanges)[subjectHits(hitsPartial)])
      newEnds <- pmin(end(queryRanges(x))[queryHits(hitsPartial)], 
                      end(searchGRanges)[subjectHits(hitsPartial)])
      for(i in 1:length(hitsPartial)){
        newAxts[[i]] <- subAln(x[queryHits(hitsPartial[i])], 
                               newStarts[i], newEnds[i],
                               select="query")
      }
    }
    # then search Axts on negative strand.
    # we need to prepare the searchGRanges on negative strand.
    searchGRanges <- GRanges(seqnames=chr,
                             ranges=IRanges(start=qSize-end+1,
                                            end=qSize-start+1),
                             strand="-")
    searchGRanges <- reduce(searchGRanges)
    hitsNegativeWithin <- findOverlaps(queryRanges(x),
                                       searchGRanges, type="within",
                                       select="all")
    indexNegativeWithin <- queryHits(hitsNegativeWithin)
    hitsNegativeAny <- findOverlaps(queryRanges(x),
                                    searchGRanges, type="any",
                                    select="all")
    indexNegativeAny <- queryHits(hitsNegativeAny)
    hitsPartial <- setdiff(hitsNegativeAny, hitsNegativeWithin)
    newAxts2 = list()
    if(length(hitsPartial) > 0L){
      newStarts <- pmax(start(queryRanges(x))[queryHits(hitsPartial)], 
                        (qSize-end(searchGRanges)+1)[subjectHits(hitsPartial)])
      newEnds <- pmin(end(queryRanges(x))[queryHits(hitsPartial)], 
                      (qSize-start(searchGRanges)+1)[subjectHits(hitsPartial)])
      for(i in 1:length(hitsPartial)){
        newAxts2[[i]] <- subAln(x[queryHits(hitsPartial[i])], 
                                newStarts[i], newEnds[i],
                                select="query")
      }
    }
    index <- sort(c(indexPositiveWithin, indexNegativeWithin))
    ans <- c(x[index], do.call(c, newAxts), do.call(c, newAxts2))
    return(ans)
  }
}



subAln <- function(x, start, end, select=c("target", "query")){
  if(length(x) != 1L){
    stop("subAln only operates on Axt of length 1")
  }
  select <- match.arg(select)
  alignedSeq1 <- strsplit(as.character(targetSeqs(x)[[1]]), "")[[1]]
  alignedSeq2 <- strsplit(as.character(querySeqs(x)[[1]]), "")[[1]]
  indexGap <- alignedSeq1 == "-" | alignedSeq1 == "." | alignedSeq1 == "_"
  seq12aln <- seq_len(length(alignedSeq1))[!indexGap]
  indexGap <- alignedSeq2 == "-" | alignedSeq2 == "." | alignedSeq1 == "_"
  seq22aln = seq_len(length(alignedSeq2))[!indexGap]
  if(select == "target"){
    alnStart <- seq12aln[start - start(targetRanges(x)) + 1L]
    alnEnd <- seq12aln[end - start(targetRanges(x)) + 1L]
  }else{
    alnStart <- seq22aln[start - start(queryRanges(x)) + 1L]
    alnEnd <- seq22aln[end - start(queryRanges(x)) + 1L]
  }
  while(alignedSeq1[alnStart] %in% c("-", ".", "_") ||
        alignedSeq2[alnStart] %in% c("-", ".", "_")){
    alnStart = alnStart + 1L
  }
  while(alignedSeq1[alnEnd] %in% c("-", ".", "_") ||
        alignedSeq2[alnEnd] %in% c("-", ".", "_")){
    alnEnd = alnEnd - 1L
  }
  ## TO DO: check newStart < newEnd
  newTargetStart <- which(alnStart == seq12aln) + start(targetRanges(x)) - 1L
  newTargetEnd <- which(alnEnd == seq12aln) + start(targetRanges(x)) - 1L
  newQueryStart <- which(alnStart == seq22aln) + start(queryRanges(x)) - 1L
  newQueryEnd <- which(alnEnd == seq22aln) + start(queryRanges(x)) - 1L
  ans <- Axt(targetRanges=GRanges(seqnames=seqnames(targetRanges(x)),
                                  ranges=IRanges(start=newTargetStart,
                                                 end=newTargetEnd),
                                  strand="+"),
             targetSeqs=DNAStringSet(paste0(alignedSeq1[alnStart:alnEnd],
                                            collapse="")),
             queryRanges=GRanges(seqnames=seqnames(queryRanges(x)),
                                 ranges=IRanges(start=newQueryStart,
                                                end=newQueryEnd),
                                 strand=strand(queryRanges(x))),
             querySeqs=DNAStringSet(paste0(alignedSeq2[alnStart:alnEnd],
                                           collapse="")),
             score=score(x), symCount=alnEnd-alnStart+1L)
  return(ans)
}


