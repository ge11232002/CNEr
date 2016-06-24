### -----------------------------------------------------------------
### match distribution plot for Axt alignment
### Expotred!
### This implementation is too slow for large Axt. Should implment in C.
setGeneric("matchDistribution", function(x, size=10000)
  standardGeneric("matchDistribution"))

setMethod("matchDistribution", signature(x="Axt"),
          function(x, size=10000){
            if(length(x) > size){
              x <- x[sort(sample.int(length(x), size))]
            }
            matchFreq <- table(unlist(mapply(paste0,
                               strsplit(as.character(targetSeqs(x)), ""), 
                               strsplit(as.character(querySeqs(x)), ""))))
            matchFreq <- matchFreq / sum(matchFreq) * 100
            freqMatrix <- matrix(0, ncol=6, nrow=6)
            colnames(freqMatrix) <- rownames(freqMatrix) <- 
              c("A", "C", "G", "T", "-", "N")
            for(i in 1:length(matchFreq)){
              splittedNames <- strsplit(names(matchFreq)[i], "")[[1]]
              freqMatrix[splittedNames[1], splittedNames[2]] <-
                matchFreq[i]
            }
            toPlot <- melt(freqMatrix)
            colnames(toPlot) <- c("target", "query", "percentage")
            ggplot(toPlot, aes_string(x="target", y="query", 
                                      fill="percentage")) +
              geom_tile() +
              theme_bw() + xlab("Target") + ylab("Query") +
              ggtitle("Distribution of matched alignments") +
              scale_fill_continuous(low="deepskyblue4", high="gold") +
              geom_text(aes(fill=toPlot$percentage, 
                            label=round(toPlot$percentage, 4)))
          }
          )

### -----------------------------------------------------------------
### summary function for Axt
### Exported!
setMethod("summary", signature=(object="Axt"),
          function(object, ...){
            x <- object
            totalLengthAlignments <- sum(as.numeric(symCount(x)))
            ## length, number of alignment,
            cat("Alignment number: ", length(x), 
                " total length: ", totalLengthAlignments, "\n")
            
            ## mismatches counts and percentage
            #compResults <- compDNAStringSet(targetSeqs(x), querySeqs(x))
            compResults <- compareStrings(targetSeqs(x), querySeqs(x))
            compResults <- table(unlist(strsplit(compResults, "")))
            #count <- sum(as.numeric(symCount(x))) - 
            #  sum(sapply(compResults, sum))
            #probability <- count / sum(as.numeric(symCount(x)))
            cat("Insertion count: ", compResults["+"],
                " percentage: ", compResults["+"] / totalLengthAlignments, "\n")
            cat("Delection count: ", compResults["-"],
                " percentage: ", compResults["-"] / totalLengthAlignments, "\n")
            cat("Mismatch count: ", compResults["?"],
                " percentage: ", compResults["?"] / totalLengthAlignments, "\n")
            cat("Match count: ", 
                sum(as.numeric(compResults[c("A", "C", "G", "T")])),
                " percentage: ", sum(as.numeric(compResults[c("A", "C", "G", "T")])) / totalLengthAlignments, "\n")
            invisible(compResults)
          }
          )

### -----------------------------------------------------------------
### dotplot for synteny from Axt object
### Exported!
dotplotAxt <- function(axt, targetSeqinfo, querySeqinfo,
                       targetChrs=NULL, queryChrs=NULL,
                       col=c("blue", "red")){
  if(!is(axt, "Axt")){
    stop("axt must be a Axt object!")
  }
  if(!is(targetSeqinfo, "Seqinfo") || !is(querySeqinfo, "Seqinfo")){
    stop("targetSeqinfo and querySeqinfo must be Seqinfo objects!")
  }
  if(!is.null(targetChrs) && !is.null(queryChrs)){
    axt <- axt[as.character(seqnames(targetRanges(axt))) %in% targetChrs & 
               as.character(seqnames(queryRanges(axt))) %in% queryChrs]
    targetSeqinfo <- targetSeqinfo[targetChrs]
    querySeqinfo <- querySeqinfo[queryChrs]
  }
  target <- targetRanges(axt)
  query <- queryRanges(axt)
  
  ## If we want to put all the segments of syntent in one plot
  ## We need to shift the coordiantes from second chromosomes in seqinfo.
  shiftCoordinatesTarget <- c(0, cumsum(as.numeric(seqlengths(targetSeqinfo)))[-length(targetSeqinfo)])
  names(shiftCoordinatesTarget) <- seqnames(targetSeqinfo)
  shiftCoordinatesQuery <- c(0, cumsum(as.numeric(seqlengths(querySeqinfo)))[-length(querySeqinfo)])
  names(shiftCoordinatesQuery) <- seqnames(querySeqinfo)
  
  startTarget <- start(target) +
    shiftCoordinatesTarget[as.character(seqnames(target))]
  ## query is kind of difficult. Because when on "-" strand,
  ## the coordinates is based on "-" strand.
  indexNegative <- as.logical(strand(query)=="-")
  startQuery <- start(query)
  startQuery[indexNegative] <- 
    seqlengths(querySeqinfo)[as.character(seqnames(query[indexNegative]))] - 
    start(query[indexNegative]) + 1
  startQuery <- startQuery + 
    shiftCoordinatesQuery[as.character(seqnames(query))]

  toPlot <- data.frame(x=startTarget, y=startQuery,
                       strand=as.factor(strand(query)))
  p <- ggplot(data=toPlot, aes_string(x="x",y="y")) +
    geom_point(aes(colour=strand), size=0.5) + theme_bw() +
    scale_colour_manual(values=col) +
    #coord_cartesian(#xlim(0, sum(as.numeric(seqlengths(targetSeqinfo)))),
    #                ylim(0, sum(as.numeric(seqlengths(querySeqinfo))))) + 
    xlab("Target") + ylab("Query") + ggtitle("Syntenic dotplot") +
    scale_x_continuous(breaks=c(0,
                                cumsum(as.numeric(seqlengths(targetSeqinfo)))),
      limits=c(0, sum(as.numeric(seqlengths(targetSeqinfo)))),
      labels=c("start", seqnames(targetSeqinfo))) +
    scale_y_continuous(breaks=c(0,
                                cumsum(as.numeric(seqlengths(querySeqinfo)))),
                       limits=c(0, sum(as.numeric(seqlengths(querySeqinfo)))),
                       labels=c("start", seqnames(querySeqinfo)))
  p
}
