### -----------------------------------------------------------------
### CNE class
### Exported!
setClass(Class="CNE",
         slots=c(assembly1="character",
                 assembly2="character",
                 window="integer",
                 identity="integer",
                 CNE12="GRangePairs",
                 CNE21="GRangePairs",
                 CNEMerged="GRangePairs",
                 CNEFinal="GRangePairs",
                 aligner="character",
                 cutoffs1="integer",
                 cutoffs2="integer",
                 smoothingWindow1="integer",
                 smoothingWindow2="integer"
         ),
         prototype=list(assembly1=character(1),
                        assembly2=character(1),
                        window=50L,
                        identity=50L,
                        CNE12=GRangePairs(),
                        CNE21=GRangePairs(),
                        CNEMerged=GRangePairs(),
                        CNEFinal=GRangePairs(),
                        aligner="blat",
                        cutoffs1=4L,
                        cutoffs2=4L,
                        smoothingWindow1=300L,
                        smoothingWindow2=300L
                        )
)

setValidity("CNE",
            function(object){
              if(length(assembly1(object)) != 1L)
                return("The filename of assembly1 must be length 1!")
              if(length(assembly2(object)) != 1L)
                return("The filename of assembly2 must be length 1!")
              if(length(object@aligner) != 1L)
                return("The aligner must be length 1!")
              if(object@identity > object@window)
                return("The identity must be equal to smaller than window")
              if(object@smoothingWindow1 > 1000 || object@smoothingWindow1 < 10)
                return("The smoothingWindow1 must be between 10 and 1000")
              if(object@smoothingWindow2 > 1000 || object@smoothingWindow2 < 10)
                return("The smoothingWindow2 must be between 10 and 1000")
              return(TRUE)
            }
)

### -----------------------------------------------------------------
### CNE class Generics
setGeneric("assembly1", function(x) standardGeneric("assembly1"))
setGeneric("assembly2", function(x) standardGeneric("assembly2"))
setGeneric("CNE12", function(x) standardGeneric("CNE12"))
setGeneric("CNE21", function(x) standardGeneric("CNE21"))
setGeneric("thresholds", function(x) standardGeneric("thresholds"))
setGeneric("CNEMerged", function(x) standardGeneric("CNEMerged"))
setGeneric("CNEFinal", function(x) standardGeneric("CNEFinal"))
setGeneric("saveCNEToSQLite", 
           function(CNE, dbName, tableName, overwrite=FALSE) 
             standardGeneric("saveCNEToSQLite"))
setGeneric("CNEDensity",
           function(dbName, tableName, assembly1, assembly2, threshold,
                    chr, start, end, windowSize, minLength=NULL)
             standardGeneric("CNEDensity"))
setGeneric("ceScan", 
           function(axts, tFilter, qFilter, qSizes, window=50L, identity=50L)
             standardGeneric("ceScan"))

### -----------------------------------------------------------------
### CNE Slot getters and setters.
### Exported!
setMethod("assembly1", "CNE", function(x) x@assembly1)
setMethod("assembly2", "CNE", function(x) x@assembly2)
setMethod("CNE12", "CNE", function(x) x@CNE12)
setMethod("CNE21", "CNE", function(x) x@CNE21)
setMethod("thresholds", "CNE", function(x)
  paste(x@identity, x@window, sep="_"))
setMethod("CNEMerged", "CNE", function(x) x@CNEMerged)
setMethod("CNEFinal", "CNE", function(x) x@CNEFinal)

### -----------------------------------------------------------------
### CNE constructor.
### Exported!
CNE <- function(assembly1=character(1), assembly2=character(1),
                window=50L, identity=50L,
                CNE12=GRangePairs(), CNE21=GRangePairs(),
                CNEMerged=GRangePairs(), CNEFinal=GRangePairs(),
                aligner="blat",
                cutoffs1=4L, cutoffs2=4L,
                smoothingWindow1=300L,
                smoothingWindow2=300L
){
  new("CNE", assembly1=assembly1, assembly2=assembly2,
      window=window, identity=identity, CNE12=CNE12, CNE21=CNE21,
      CNEMerged=CNEMerged, CNEFinal=CNEFinal,
      aligner=aligner, cutoffs1=cutoffs1, cutoffs2=cutoffs2,
      smoothingWindow1=smoothingWindow1, 
      smoothingWindow2=smoothingWindow2)
}