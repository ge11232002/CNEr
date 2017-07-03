### -----------------------------------------------------------------
### Process the direct output from goseq package.
### It outputs the data.frame with indent "."
### -----------------------------------------------------------------
getChildTerms <- function(x, subset, goRelatives, indent="", childEnvir){
  result <- character()
  if (is.element(x, subset)){
    term <- getGOTerm(x)[[1]]
    displayLabel <- paste(indent, substr(term, 1, 60), sep="")
    result[displayLabel] <- x
    ## this is the modification to make the table non-redundant!!
    subset <- setdiff(subset, x) 
  }
  kids <- intersect(get(x, envir=childEnvir), goRelatives)
  indent <- paste(indent, ".")
  for (kid in kids){
    kidResult <- getChildTerms(kid, subset, goRelatives,
                               indent=indent, childEnvir)
    subset <- setdiff(subset, kidResult) 
    ## again the non-redundancy modification
    result <- c(result, kidResult)
  }
  result
}

### -----------------------------------------------------------------
### make the nice table of GO
### Need to library(GO.db) first
### Not exported!
hierarchicGOTable <- function(x, onto=c("BP", "CC", "MF"),
                              maxNumberOfTerms=100){
  onto <- match.arg(onto)
  # x is the output from the goseq pipeline
  rownames(x) <- x$category
  x <- x[x$ontology==onto, c("over_represented_pvalue", "numDEInCat", 
                             "numInCat", "term", "category")]
  colnames(x)[1:3] <- c("Pvalue", "Count", "Size")
  if (!is.data.frame(x)){
    return(data.frame())
  }
  if (nrow(x) > maxNumberOfTerms){
    x = x[1:maxNumberOfTerms, ]
  }
  if (nrow(x) == 0){
    return(data.frame())
  }
  
  if (onto == "CC"){
    ANCESTOR <- GOCCANCESTOR
    OFFSPRING <- GOCCOFFSPRING
    CHILDREN <- GOCCCHILDREN
  }
  if (onto == "BP"){
    ANCESTOR <- GOBPANCESTOR
    OFFSPRING <- GOBPOFFSPRING
    CHILDREN <- GOBPCHILDREN
  }
  if (onto == "MF"){
    ANCESTOR <- GOMFANCESTOR
    OFFSPRING <- GOMFOFFSPRING
    CHILDREN <- GOMFCHILDREN
  }
  
  goIds <- rownames(x)
  goAncestorList <- mget(goIds, envir=ANCESTOR)
  goRoots <- character()
  for (goId in goIds){
    if (length(intersect(goIds, goAncestorList[[goId]])) == 0){
      goRoots[goId] <- goId
    }
  }
  goOffsprings <- unique(unlist(mget(goIds, envir=OFFSPRING)))
  goAncestors <- unique(unlist(goAncestorList))
  goRelatives <- union(intersect(goAncestors, goOffsprings), goIds)
  
  tables <- list()
  for (i in 1:length(goRoots)){
    childTerms <- getChildTerms(goRoots[i], goIds, goRelatives, indent="",
                                CHILDREN)
    tables[[i]] <-
      data.frame(terms=sub("^ ", "", names(childTerms)),
                 pValues=format(x[childTerms, "Pvalue"], digits=3,
                                scientific=TRUE),
                 counts=paste(x[childTerms, "Count"], x[childTerms, "Size"],
                              sep="/"),
                 category=x[childTerms, "category"],
                 stringsAsFactors=FALSE
      )
  }
  ans <- do.call(rbind, tables)
  return(ans)
}

## ------------------------------------------------------------------
## Find the closest root of GO terms
## ------------------------------------------------------------------
findGORoot <- function(goIds, onto=c("BP", "CC", "MF")){
  onto <- match.arg(onto)
  
  if (onto == "CC"){
    ANCESTOR = GOCCANCESTOR
    OFFSPRING = GOCCOFFSPRING
    CHILDREN = GOCCCHILDREN
  }
  if (onto == "BP"){
    ANCESTOR = GOBPANCESTOR
    OFFSPRING = GOBPOFFSPRING
    CHILDREN = GOBPCHILDREN
  }
  if (onto == "MF"){
    ANCESTOR = GOMFANCESTOR
    OFFSPRING = GOMFOFFSPRING
    CHILDREN = GOMFCHILDREN
  }
  goAncestorList = mget(goIds, envir=ANCESTOR)
  goTable <- table(unlist(goAncestorList))
  names(goTable)[goTable == length(goAncestorList)]
}

### -----------------------------------------------------------------
### Read the GAF file
### -----------------------------------------------------------------
readGAF <- function(fn){
  lines <- readLines(fn)
  lines <- lines[-(1:3)]
  geneIDs <- sapply(strsplit(lines, "\t"), "[", 2)
  GOTerms <- sapply(strsplit(lines, "\t"), "[", 5)
  ans <- split(GOTerms, geneIDs)
  return(ans)
}

### -----------------------------------------------------------------
### Add ancestor GO
### Exported!
addAncestorGO <- function(go){
  if(!is(go, "list")){
    stop("`go` must be a list!")
  }
  goID2Ancestor <- c(as.list(GOBPANCESTOR), as.list(GOMFANCESTOR), 
                     as.list(GOCCANCESTOR))
  allGo <- unlist(go)
  goID2Ancestor <- lapply(goID2Ancestor, function(x){x[x%in%allGo]})
  newGo <- lapply(relist(mapply(append, allGo, goID2Ancestor[allGo]), go),
                  unlist)
  newGo <- lapply(newGo, function(x){if(is.null(x)){character(0)}else{x}})
  return(newGo)
}
