### -----------------------------------------------------------------
### orgKEGGIds2EntrezIDs: This script is supposed to parse the html page of 
###              certain species's pathway from KEGG and 
###              download the associated information for each KEGG pathway ID.
### Exported!
orgKEGGIds2EntrezIDs <- function(organism="Homo sapiens"){
  ## Species mapping
  organismMapping <- read.table("http://rest.kegg.jp/list/organism",
                                header=FALSE, sep="\t", quote="",
                                comment.char="")
  organismID <- organismMapping[grepl(organism, organismMapping$V3,
                                      ignore.case=TRUE),
                                2, drop=TRUE]
  if(length(organismID) == 0L){
    stop("The provided organism is not available.",
         "Please refer to http://rest.kegg.jp/list/organism for available",
         "organisms")
  }
  html <- readLines(paste0("http://www.genome.jp/kegg-bin/show_organism?menu_type=pathway_maps&org=", organismID))
  html <- grep("^\\d{5}&", html, value=TRUE)
  ### Hopefully the ID is always 5-digit
  pathwayIDs <- paste0(organismID, substr(html, 1L, 5L))
  groups <- sample(rep_len(1L:ceiling(length(pathwayIDs) / 10),
                           length.out=length(pathwayIDs)))
  pathwayIDs <- split(pathwayIDs, groups)
  
  ## query with KEGG Rest server with 10 entries (maximal) a time, 
  ## and 200s to 400s gap between each query.
  query <- lapply(pathwayIDs,
                  function(x){Sys.sleep(sample(200L:400L, size=1L));keggGet(x)})
  
  ## re-organise the query object
  pathways <- list()
  for(i in 1:length(query)){
    for(j in 1:length(query[[i]])){
      pathways[[query[[i]][[j]]$ENTRY]] <-
        query[[i]][[j]]
    }
  }
  ## Get the Pathway IDs to Entrez Gene IDs mapping
  pathwayIDs2GeneIDs <- list()
  for(i in 1:length(pathways)){
    genesInfo <- pathways[[i]]$GENE
    if(is.null(genesInfo)){
      pathwayIDs2GeneIDs[[pathways[[i]]$ENTRY]] <- NA
      next
    }
    pathwayIDs2GeneIDs[[pathways[[i]]$ENTRY]] <-
      genesInfo[seq(1, length(genesInfo), by=2)]
  }
  pathwayIDs2GeneIDs <- pathwayIDs2GeneIDs[!is.na(pathwayIDs2GeneIDs)]
  return(pathwayIDs2GeneIDs)
}