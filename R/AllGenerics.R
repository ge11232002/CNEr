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
### general
### Not Exported!
setGeneric("clone", function(x, ...) standardGeneric("clone"))

### -----------------------------------------------------------------
### ceScan
### Exported!
setGeneric("ceScan", 
           function(axts, tFilter, qFilter, qSizes, thresholds="49_50")
           standardGeneric("ceScan"))


