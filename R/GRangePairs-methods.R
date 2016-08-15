setGeneric("syntenicDotplot", 
           function(x, firstSeqlengths=NULL, secondSeqlengths=NULL,
                    firstChrs=NULL, secondChrs=NULL,
                    col=c("blue", "red"), type=c("line", "dot"))
  standardGeneric("syntenicDotplot"))

setMethod("syntenicDotplot", signature=(x="GRangePairs"),
          function(x, firstSeqlengths=NULL, secondSeqlengths=NULL,
                   firstChrs=NULL, secondChrs=NULL,
                   col=c("blue", "red"), type=c("line", "dot")){
            syntenicPlotGRangePairs(x, firstSeqlengths=firstSeqlengths, 
                                    secondSeqlengths=secondSeqlengths,
                                    firstChrs=firstChrs,
                                    secondChrs=secondChrs,
                                    col=col, type=type)
          })

syntenicPlotGRangePairs <- function(x, firstSeqlengths=NULL, 
                                    secondSeqlengths=NULL,
                                    firstChrs=NULL, secondChrs=NULL,
                                    col=c("blue", "red"),
                                    type=c("line", "dot")){
  type <- match.arg(type)
  
  if(is.null(firstSeqlengths)){
    targetSeqlengths <- seqlengths(first(x))
    if(any(is.na(targetSeqlengths))){
      stop("When firstSeqlengths is NULL, the seqlengths must exist in x!")
    }
  }else{
    targetSeqlengths <- firstSeqlengths
  }
  
  if(is.null(secondSeqlengths)){
    querySeqlengths <- seqlengths(second(x))
    if(any(is.na(querySeqlengths))){
      stop("When secondSeqlengths is NULL, the seqlengths must exist in x!")
    }
  }else{
    querySeqlengths <- secondSeqlengths
  }
  
  if(!is.null(firstChrs)){
    x <- x[seqnames(first(x)) %in% firstChrs]
    targetSeqlengths <- targetSeqlengths[firstChrs]
  }
  if(!is.null(secondChrs)){
    x <- x[seqnames(second(x)) %in% secondChrs]
    querySeqlengths <- querySeqlengths[secondChrs]
  }
  target <- first(x)
  query <- second(x)
  
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
  if(type == "line"){
    p <- ggplot(data=toPlot, aes_string(x="x",y="y", xend="xEnd", yend="yEnd")) +
      geom_segment(aes(colour=strand)) + theme_bw() +
      scale_colour_manual(values=col) +
      xlab("First") + ylab("Second") + ggtitle("Syntenic dotplot") +
      scale_x_continuous(breaks=c(0,
                                  cumsum(as.numeric(targetSeqlengths))),
                         limits=c(0, sum(as.numeric(targetSeqlengths))),
                         labels=c("start", names(targetSeqlengths))) +
      scale_y_continuous(breaks=c(0,
                                  cumsum(as.numeric(querySeqlengths))),
                         limits=c(0, sum(as.numeric(querySeqlengths))),
                         labels=c("start", names(querySeqlengths)))
  }else if(type == "dot"){
    p <- ggplot(data=toPlot, aes_string(x="x",y="y")) +
      geom_point(aes(colour=strand), size=0.5) + theme_bw() +
      scale_colour_manual(values=col) +
      xlab("First") + ylab("Second") + ggtitle("Syntenic dotplot") +
      scale_x_continuous(breaks=c(0,
                                  cumsum(as.numeric(targetSeqlengths))),
                         limits=c(0, sum(as.numeric(targetSeqlengths))),
                         labels=c("start", names(targetSeqlengths))) +
      scale_y_continuous(breaks=c(0,
                                  cumsum(as.numeric(querySeqlengths))),
                         limits=c(0, sum(as.numeric(querySeqlengths))),
                         labels=c("start", names(querySeqlengths)))
  }
  p
}