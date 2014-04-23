### -----------------------------------------------------------------
### summary function for Axt
### Not Expotred!
### This implementation is too slow for large Axt. Should implment in C.
setMethod("matchDistr", signature(x="Axt"),
          function(x){
            matchFreq <- table(unlist(mapply(paste0, 
                               strsplit(as.character(targetSeqs(x)), ""), 
                               strsplit(as.character(querySeqs(x)), ""))))
            matchFreq <- matchFreq / sum(matchFreq)
            freqMatrix <- matrix(0, ncol=6, nrow=6)
            colnames(freqMatrix) <- c("A", "C", "G", "T", "-", "N")
            rownames(freqMatrix) <- c("A", "C", "G", "T", "-", "N")
            for(i in 1:length(matchFreq)){
              splittedNames <- strsplit(names(matchFreq)[i], "")[[1]]
              freqMatrix[splittedNames[1], splittedNames[2]] <-
                matchFreq[i]
            }
            return(freqMatrix)
          }
          )

