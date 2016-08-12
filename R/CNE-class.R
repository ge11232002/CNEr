### -----------------------------------------------------------------
### CNE class
### Exported!
setClass(Class="CNE",
         slots=c(assembly1Fn="character",
                 assembly2Fn="character",
                 axt12Fn="character",
                 axt21Fn="character",
                 window="integer",
                 identity="integer",
                 CNE12="GRangePairs",
                 CNE21="GRangePairs",
                 CNEMerged="GRangePairs",
                 CNEFinal="GRangePairs",
                 aligner="character",
                 cutoffs1="integer",
                 cutoffs2="integer"
         ),
         prototype=list(assembly1Fn=character(1),
                        assembly2Fn=character(1),
                        axt12Fn=character(1),
                        axt21Fn=character(1),
                        window=50L,
                        identity=50L,
                        CNE12=GRangePairs(),
                        CNE21=GRangePairs(),
                        CNEMerged=GRangePairs(),
                        CNEFinal=GRangePairs(),
                        aligner="blat",
                        cutoffs1=4L,
                        cutoffs2=4L
                        )
)

setValidity("CNE",
            function(object){
              if(length(object@assembly1Fn) != 1L)
                return("The filename of assembly1Fn must be length 1!")
              if(length(object@assembly2Fn) != 1L)
                return("The filename of assembly2Fn must be length 1!")
              if(length(object@aligner) != 1L)
                return("The aligner must be length 1!")
              if(object@identity > object@window)
                return("The identity must be equal to smaller than window")
              return(TRUE)
            }
)

### -----------------------------------------------------------------
### CNE class Generics
#setGeneric("assembly1", function(x) standardGeneric("assembly1"))
#setGeneric("assembly2", function(x) standardGeneric("assembly2"))
setGeneric("CNE12", function(x) standardGeneric("CNE12"))
setGeneric("CNE21", function(x) standardGeneric("CNE21"))
setGeneric("thresholds", function(x) standardGeneric("thresholds"))
setGeneric("CNEMerged", function(x) standardGeneric("CNEMerged"))
setGeneric("CNEFinal", function(x) standardGeneric("CNEFinal"))
#setGeneric("saveCNEToSQLite", 
#           function(CNE, dbName, tableName, overwrite=FALSE) 
#             standardGeneric("saveCNEToSQLite"))
setGeneric("ceScan", 
           function(axts, tFilter, qFilter, qSizes, window=50L, identity=50L)
             standardGeneric("ceScan"))

### -----------------------------------------------------------------
### CNE Slot getters and setters.
### Exported!
#setMethod("assembly1", "CNE", function(x) x@assembly1)
#setMethod("assembly2", "CNE", function(x) x@assembly2)
setMethod("CNE12", "CNE", function(x) x@CNE12)
setMethod("CNE21", "CNE", function(x) x@CNE21)
setMethod("thresholds", "CNE", function(x)
  paste(x@identity, x@window, sep="_"))
setMethod("CNEMerged", "CNE", function(x) x@CNEMerged)
setMethod("CNEFinal", "CNE", function(x) x@CNEFinal)

### -----------------------------------------------------------------
### CNE constructor.
### Exported!
CNE <- function(assembly1Fn=character(1), assembly2Fn=character(1),
                axt12Fn=character(1), axt21Fn=character(1),
                window=50L, identity=50L,
                CNE12=GRangePairs(), CNE21=GRangePairs(),
                CNEMerged=GRangePairs(), CNEFinal=GRangePairs(),
                aligner="blat",
                cutoffs1=4L, cutoffs2=4L
){
  new("CNE", assembly1Fn=assembly1Fn, assembly2Fn=assembly2Fn,
      axt12Fn=axt12Fn, axt21Fn=axt21Fn,
      window=window, identity=identity, CNE12=CNE12, CNE21=CNE21,
      CNEMerged=CNEMerged, CNEFinal=CNEFinal,
      aligner=aligner, cutoffs1=cutoffs1, cutoffs2=cutoffs2)
}
