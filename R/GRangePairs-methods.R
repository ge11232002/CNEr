setGeneric("syntenicDotplot", 
           function(x, firstSeqlengths=NULL, lastSeqlengths=NULL,
                    firstChrs=NULL, lastChrs=NULL,
                    col=c("blue", "red"))
  standardGeneric("syntenicDotplot"))

setMethod("syntenicDotplot", signature=(x="GRangePairs"),
          function(x, firstSeqlengths=NULL, lastSeqlengths=NULL,
                   firstChrs=NULL, lastChrs=NULL,
                   col=c("blue", "red")){
            syntenicPlotGRangePairs(x, firstSeqlengths=firstSeqlengths, 
                                    lastSeqlengths=lastSeqlengths,
                                    firstChrs=firstChrs,
                                    lastChrs=lastChrs,
                                    col=col)
          })

syntenicPlotGRangePairs <- function(x, firstSeqlengths=NULL, 
                                    lastSeqlengths=NULL,
                                    firstChrs=NULL, lastChrs=NULL,
                                    col=c("blue", "red")){
  if(is.null(firstSeqlengths)){
    targetSeqlengths <- seqlengths(first(x))
  }
  if(is.null(lastSeqlengths)){
    querySeqlengths <- seqlengths(last(x))
  }
  target <- first(x)
  query <- last(x)
  
  if(!is.null(firstChrs)){
    x <- x[as.character(seqnames(first(x))) %in% firstChrs]
    targetSeqlengths <- targetSeqlengths[firstChrs]
  }
  if(!is.null(lastChrs)){
    x <- x[as.character(seqnames(last(x))) %in% lastChrs]
    querySeqlengths <- querySeqlengths[lastChrs]
  }
  target <- first(x)
  query <- last(x)
  
  ## If we want to put all the segments of syntent in one plot
  ## We need to shift the coordiantes from second chromosomes in seqlengths
  shiftCoordinatesTarget <- c(0, cumsum(as.numeric(targetSeqlengths))[-length(targetSeqlengths)])
  names(shiftCoordinatesTarget) <- names(targetSeqlengths)
  shiftCoordinatesQuery <- c(0, cumsum(as.numeric(querySeqlengths))[-length(querySeqlengths)])
  names(shiftCoordinatesQuery) <- names(querySeqlengths)
  
  startTarget <- start(target) +
    shiftCoordinatesTarget[as.character(seqnames(target))]
  endTarget <- end(target) + 
    shiftCoordinatesTarget[as.character(seqnames(target))]
  startQuery <- start(query) +
    shiftCoordinatesQuery[as.character(seqnames(query))]
  endQuery <- end(query) +
    shiftCoordinatesQuery[as.character(seqnames(query))]
  combinedStrands <- rep("+", length(target))
  combinedStrands[xor(as.character(strand(target))=="+", 
                      as.character(strand(query))=="+")] <- "-"
  
  toPlot <- data.frame(x=startTarget, xEnd=endTarget,
                       y=startQuery, yEnd=endQuery,
                       strand=factor(combinedStrands, levels=c("+", "-")))
  p <- ggplot(data=toPlot, aes_string(x="x",y="y", xend="xEnd", yend="yEnd")) +
    geom_segment(aes(colour=strand)) + theme_bw() +
    scale_colour_manual(values=col) +
    xlab("First") + ylab("Last") + ggtitle("Syntenic dotplot") +
    scale_x_continuous(breaks=c(0,
                                cumsum(as.numeric(targetSeqlengths))),
                       limits=c(0, sum(as.numeric(targetSeqlengths))),
                       labels=c("start", names(targetSeqlengths))) +
    scale_y_continuous(breaks=c(0,
                                cumsum(as.numeric(querySeqlengths))),
                       limits=c(0, sum(as.numeric(querySeqlengths))),
                       labels=c("start", names(querySeqlengths)))
  p
  
}