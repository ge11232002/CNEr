### -----------------------------------------------------------------
### match distribution plot for Axt alignment
### Expotred!
### This implementation is too slow for large Axt. Should implment in C.
setGeneric("matchDistribution", function(x, size=10000, title=NULL)
  standardGeneric("matchDistribution"))

setMethod("matchDistribution", signature(x="Axt"),
          function(x, size=10000, title=NULL){
            axtMatchDistribution(x, size=size, title=title)
          }
          )

axtMatchDistribution <- function(x, size=10000, title=NULL){
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
    ## Some ambiguous bases are not considered here.
    if(!all(splittedNames %in% colnames(freqMatrix))){
      next
    }
    freqMatrix[splittedNames[1], splittedNames[2]] <-
      matchFreq[i]
  }
  toPlot <- melt(freqMatrix)
  colnames(toPlot) <- c("target", "query", "percentage")
  ggplot(toPlot, aes_string(x="target", y="query", 
                            fill="percentage")) +
    geom_tile() +
    theme_bw() + xlab("Target") + ylab("Query") +
    ggtitle(ifelse(is.null(title), "Distribution of matched bases",
                  title)) +
    scale_fill_continuous(low="deepskyblue4", high="gold") +
    geom_text(aes(label=round(toPlot$percentage, 1)))
}

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
setMethod("syntenicDotplot", signature=(x="Axt"),
          function(x, firstSeqlengths=NULL, secondSeqlengths=NULL,
                   firstChrs=NULL, secondChrs=NULL,
                   col=c("blue", "red"),
                   type=c("line", "dot")){
            if(!is.null(firstSeqlengths)){
              seqlengths(first(x)) <- 
                firstSeqlengths[names(seqlengths(first(x)))]
            }
            if(!is.null(secondSeqlengths)){
              seqlengths(second(x)) <- 
                secondSeqlengths[names(seqlengths(second(x)))]
            }
            axt <- fixCoordinates(x)
            class(axt) <- "GRangePairs"
            syntenicPlotGRangePairs(axt, firstSeqlengths=NULL,
                                    secondSeqlengths=NULL,
                                    firstChrs=firstChrs,
                                    secondChrs=secondChrs,
                                    col=col, type=type)
          })

### -----------------------------------------------------------------
### makeAxtTracks
### Exported!!
makeAxtTracks <- function(x){
  if(!is(x, "Axt")){
    stop(deparse(substitute(x)), " must be a `Axt`` class!")
  }
  x <- fixCoordinates(x)
  
  targetAxt <- targetRanges(x)
  queryAxt <- queryRanges(x)
  
  targetAxt$name <- as.character(queryAxt)
  queryAxt$name <- as.character(targetAxt)
  
  export.bed(targetAxt, "targetAxt.bed")
  export.bed(queryAxt, "queryAxt.bed")
  
  invisible(list(targetAxt, queryAxt))
}
