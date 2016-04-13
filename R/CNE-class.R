### -----------------------------------------------------------------
### CNE class
### Exported!
setClass(Class="CNE",
         slots=c(assembly1="character",
                 assembly2="character",
                 thresholds="character",
                 CNE1="list",
                 CNE2="list",
                 CNEMerged="list",
                 CNERepeatsFiltered="list",
                 alignMethod="character"
         )
)

setValidity("CNE",
            function(object){
              if(length(object@assembly1) != 1L)
                return("The name of assembly1 must be length 1!")
              if(length(object@alignMethod) != 1L)
                return("The align method must be length 1!")
              if(length(object@assembly2) != 1L)
                return("The name of assembly2 must be length 1!")
              if(!all(grepl("^\\d+_\\d+$", object@thresholds)))
                return("The thresholds must be in format of 49_50!")
              if(any(as.integer(
                sapply(strsplit(object@thresholds, "_"), "[", 2))
                < as.integer(
                  sapply(strsplit(object@thresholds, "_"), "[", 1))))
                return("The window size cannot be smaller than identity score!")
              if(length(object@CNE1) != length(object@thresholds) ||
                 length(object@CNE2) != length(object@thresholds) ||
                 length(object@CNEMerged) != length(object@thresholds) ||
                 length(object@CNERepeatsFiltered) != length(object@thresholds))
                return("The number of cne tables must be same with
                       number of thresholds!")
              return(TRUE)
            }
)

### -----------------------------------------------------------------
### CNE class related
###
setGeneric("assembly1", function(x) standardGeneric("assembly1"))
setGeneric("assembly2", function(x) standardGeneric("assembly2"))
setGeneric("CNE1", function(x) standardGeneric("CNE1"))
setGeneric("CNE2", function(x) standardGeneric("CNE2"))
setGeneric("thresholds", function(x) standardGeneric("thresholds"))
setGeneric("CNEMerged", function(x) standardGeneric("CNEMerged"))
setGeneric("CNERepeatsFiltered", function(x) 
  standardGeneric("CNERepeatsFiltered"))

setGeneric("saveCNEToSQLite", 
           function(CNE, dbName, tableName, overwrite=FALSE) 
             standardGeneric("saveCNEToSQLite"))
setGeneric("CNEDensity",
           function(dbName, tableName, assembly1, assembly2, threshold,
                    chr, start, end, windowSize, minLength=NULL)
             standardGeneric("CNEDensity"))

### -----------------------------------------------------------------
### ceScan
### Exported!
setGeneric("ceScan", 
           function(axts, tFilter, qFilter, qSizes, thresholds="49_50")
             standardGeneric("ceScan"))

### -----------------------------------------------------------------
### CNE Slot getters and setters.
### Exported!
setMethod("assembly1", "CNE", function(x) x@assembly1)
setMethod("assembly2", "CNE", function(x) x@assembly2)
setMethod("CNE1", "CNE", function(x) x@CNE1)
setMethod("CNE2", "CNE", function(x) x@CNE2)
setMethod("thresholds", "CNE", function(x) x@thresholds)
setMethod("CNEMerged", "CNE", function(x) x@CNEMerged)
setMethod("CNERepeatsFiltered", "CNE", function(x) x@CNERepeatsFiltered)

### -----------------------------------------------------------------
### CNE constructor.
### Exported!
CNE <- function(assembly1=character(), assembly2=character(),
                thresholds=character(),
                CNE1=list(), CNE2=list(),
                CNEMerged=list(), CNERepeatsFiltered=list(),
                alignMethod=character()
){
  new("CNE", assembly1=assembly1, assembly2=assembly2,
      thresholds=thresholds, CNE1=CNE1, CNE2=CNE2,
      CNEMerged=CNEMerged, CNERepeatsFiltered=CNERepeatsFiltered,
      alignMethod=alignMethod)
}