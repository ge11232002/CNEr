### -----------------------------------------------------------------
### chromosome to colour mapping
### Not exported!
chr2colour <- function(chromosomes){
  ### diverging colour scheme from http://colorbrewer2.org
  divergingColour <- c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b",
                       "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd",
                       "#5e4fa2")
  uniqChromosomes <- as.character(unique(chromosomes))
  chrColours <- rep(divergingColour,
                    ceiling(length(uniqChromosomes) / length(divergingColour)))
  mapping <- chrColours[1:length(uniqChromosomes)]
  names(mapping) <- uniqChromosomes
  return(mapping[as.character(chromosomes)])
}

### -----------------------------------------------------------------
### readAncoraCNEs: this is for internal use in Lenhard group.
### It reads the CNE file (filename format: cne2wBf_hg38_mm10_50_50)
### and return a GRanges.
### Exported!
readAncora <- function(fn, assembly=NULL,
                       tAssemblyFn=NULL, qAssemblyFn=NULL){
  assembly1 <- strsplit(basename(fn), split="_")[[1]][2]
  assembly2 <- strsplit(basename(fn), split="_")[[1]][3]
  cne <- read.delim(fn, header=FALSE)
  
  # Prepare the seqinfo when available
  seqinfoTarget <- NULL
  if(!is.null(tAssemblyFn)){
    seqinfoTarget <- seqinfoFn(tAssemblyFn)
  }
  seqinfoQuery <- NULL
  if(!is.null(qAssemblyFn)){
    seqinfoQuery <- seqinfoFn(qAssemblyFn)
  }
  
  ans <- GRangePairs(first=GRanges(seqnames=cne[[1]],
                                   ranges=IRanges(start=cne[[2]]+1,
                                                  end=cne[[3]]),
                                   strand="*",
                                   name=paste0(cne[[4]], ":", 
                                               (cne[[5]]+1), "-", cne[[6]]),
                                   itemRgb=chr2colour(cne[[4]]),
                                   seqinfo=seqinfoTarget),
                     second=GRanges(seqnames=cne[[4]],
                                  ranges=IRanges(start=cne[[5]]+1,
                                                 end=cne[[6]]),
                                  strand="*",
                                  name=paste0(cne[[1]], ":", 
                                              (cne[[2]]+1), "-", cne[[3]]),
                                  itemRgb=chr2colour(cne[[1]]),
                                  seqinfo=seqinfoQuery)
                     )
  if(is.null(assembly)){
    ## Real both assemblies
    return(ans)
  }else if(assembly == assembly1){
    ## Read first assembly
    return(first(ans))
  }else if(assembly == assembly2){
    ## Read last assembly
    return(second(ans))
  }else{
    stop("Wrongly specified assembly!")
  }
}

### -----------------------------------------------------------------
### readAncoraIntoSQLite: read the Ancora format CNE files into a SQLite
### Exported!
readAncoraIntoSQLite <- function(cneFns, dbName, overwrite=FALSE){
  tableNames <- c()
  for(cneFn in cneFns){
    ## cneFn in format of cne2wBf_AstMex102_danRer10_48_50
    tableName <- sub("^cne2wBf_", "", basename(cneFn))
    tableNames <- c(tableNames, tableName)
    df <- readAncora(cneFn, assembly=NULL)
    saveCNEToSQLite(df, dbName=dbName, tableName=tableName,
                    overwrite=overwrite)
  }
  invisible(tableNames)
}

### -----------------------------------------------------------------
### makeCNEDensity: make the Ancora downloads-like bed files, bedgraph, 
### bigwig files
###   from a GrangePairs of CNEs.
### Exported!
makeCNEDensity <- function(x, outputDir=".",
                           genomeFirst="first", genomeSecond="second",
                           threshold="50_50",
                           windowSizeFirst=300, ## kb
                           windowSizeSecond=300 ## kb
                          ){
  if(!is(x, "GRangePairs")){
    stop("`x` must be a GRangePairs object!")
  }
  if(seqlengthsNA(x)){
    stop("seqlengths must be provided in `x`!")
  }
  if(length(x) == 0L){
    warning("No CNEs in `x`!")
    return(FALSE)
  }
  ## make the bed files
  message("Making bed files...")
  bedFirst <- first(x)
  bedSecond <- second(x)
  mcols(bedFirst) <- DataFrame(name=as.character(bedSecond),
                               score=0,
                               itemRgb=chr2colour(as.character(
                                 seqnames(bedSecond)))
                               )
  mcols(bedSecond) <- DataFrame(name=as.character(bedFirst),
                                score=0,
                                itemRgb=chr2colour(as.character(
                                  seqnames(bedFirst)))
                                )
  firstTrackLine <- new("BasicTrackLine", itemRgb=TRUE,
                        name=paste(genomeFirst, "CNEs", threshold),
                        description=paste(genomeFirst, "CNEs", threshold)
                        )
  secondTrackLine <- new("BasicTrackLine", itemRgb=TRUE,
                         name=paste(genomeSecond, "CNEs", threshold),
                         description=paste(genomeSecond, "CNEs", threshold)
                         )
  bedFnFirst <- file.path(outputDir, 
                          paste0("CNE_", genomeFirst, "_",
                                 genomeSecond, "_",
                                 threshold, ".bed"))
  export.bed(bedFirst, con=bedFnFirst,
             trackLine=firstTrackLine)
  bedFnSecond <- file.path(outputDir, 
                           paste0("CNE_", genomeSecond, "_",
                                  genomeFirst, "_",
                                  threshold, ".bed"))
  
  export.bed(bedSecond, con=bedFnSecond,
             trackLine=secondTrackLine)
  
  # Make the bedGraph files
  message("Making bedGraph files...")
  bedFirst <- reduce(bedFirst, ignore.strand=TRUE)
  covFirst <- coverage(bedFirst)
  densityFirst <- runmean(covFirst, k=windowSizeFirst*1000, 
                          endrule = "constant") * 100
  
  bedSecond <- reduce(bedSecond, ignore.strand=TRUE)
  covSecond <- coverage(bedSecond)
  densitySecond <- runmean(covSecond, k=windowSizeSecond*1000,
                           endrule = "constant") * 100
 
  firstTrackLine <- new("GraphTrackLine",
                        name=paste(genomeFirst, "CNEs density", threshold),
                        description=paste(genomeFirst, "CNEs density", 
                                          threshold),
                        visibility="full", type="bedGraph", autoScale=TRUE
                        )
  secondTrackLine <- new("GraphTrackLine",
                        name=paste(genomeSecond, "CNEs density", threshold),
                        description=paste(genomeSecond, "CNEs density", 
                                          threshold),
                        visibility="full", type="bedGraph", autoScale=TRUE
  )
  bwFnFirst <- file.path(outputDir, 
                         paste0("CNE_density_", genomeFirst, "_",
                                genomeSecond, "_",
                                threshold, ".bedGraph"))
  export.bedGraph(densityFirst, con=bwFnFirst, trackLine=firstTrackLine)
  
  bwFnSecond <- file.path(outputDir,
                          paste0("CNE_density_", genomeSecond, "_",
                                 genomeFirst, "_",
                                 threshold, ".bedGraph"))
  export.bedGraph(densitySecond, con=bwFnSecond, trackLine=secondTrackLine)
  
  # Make bigwig files
  message("Making bigwig files...")
  bigwigFnFirst <- file.path(outputDir, 
                             paste0("CNE_density_", genomeFirst, "_",
                                    genomeSecond, "_", threshold, ".bw"))
  export.bw(densityFirst, con=bwFnFirst)
  bigwigFnSecond <- file.path(outputDir, 
                              paste0("CNE_density_", genomeSecond, "_",
                                     genomeFirst, "_", threshold, ".bw"))
  export.bw(densitySecond, con=bigwigFnSecond)
  
  invisible(c(bedFnFirst, bedFnSecond, bwFnFirst, bwFnSecond,
              bigwigFnFirst, bigwigFnSecond))
}

### -----------------------------------------------------------------
### makeAncoraFiles: in the format of cne2wBf_GmorY1_dm6_21_30 for loading
###   into Ancora.
### Exported!!!
makeAncoraFiles <- function(cne, outputDir=".",
                            genomeFirst="first", genomeSecond="second",
                            threshold="50_50"){
  # cne is a GRangePairs object
  Sys.setlocale("LC_COLLATE","C")
  unsortedNames <- c(genomeFirst, genomeSecond)
  sortedNames <- sort(c(genomeFirst, genomeSecond))
  fileName <- paste("cne2wBf", paste(sortedNames, collapse="_"), 
                    threshold, sep="_")
  if(identical(sortedNames, unsortedNames)){
    ## The order is right
    ans <- data.frame(seqnames(first(cne)),
                      start(first(cne))-1L, end(first(cne)),
                      seqnames(second(cne)),
                      start(second(cne))-1L, end(second(cne)),
                      strand(second(cne)),
                      mcols(cne)$score,
                      mcols(cne)$cigar)
  }else{
    ## The order is not  
    ans <- data.frame(seqnames(second(cne)),
                      start(second(cne))-1L, end(second(cne)),
                      seqnames(first(cne)),
                      start(first(cne))-1L, end(first(cne)),
                      strand(first(cne)),
                      mcols(cne)$score,
                      mcols(cne)$cigar)
  }
  dir.create(outputDir, showWarnings=FALSE, recursive=TRUE)
  write.table(ans, file=file.path(outputDir, fileName), 
              sep="\t", quote=FALSE, row.names=FALSE,
              col.names=FALSE)
  invisible(file.path(outputDir, fileName))
}