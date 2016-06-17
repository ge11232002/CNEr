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
setMethod("mismatchSummary", signature(x="Axt"),
          function(x, ...){
            compResults <- compDNAStringSet(targetSeqs(x), querySeqs(x))
            count <- sum(as.numeric(symCount(x))) - 
                       sum(sapply(compResults, sum))
            probability <- count / sum(as.numeric(symCount(x)))
            return(c("Count"=count, "Probability"=probability))
          }
          )