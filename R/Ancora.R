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
readAncora <- function(fn, assembly=NULL){
  assembly1 <- strsplit(basename(fn), split="_")[[1]][2]
  assembly2 <- strsplit(basename(fn), split="_")[[1]][3]
  cne <- read_tsv(fn, col_names=FALSE)
  ans <- GRangePairs(first=GRanges(seqnames=cne$X1,
                                   ranges=IRanges(start=cne$X2+1,
                                                  end=cne$X3),
                                   strand="*",
                                   name=paste0(cne$X4, ":", 
                                               (cne$X5+1), "-", cne$X6),
                                   itemRgb=chr2colour(cne$X4)),
                     second=GRanges(seqnames=cne$X4,
                                  ranges=IRanges(start=cne$X5+1,
                                                 end=cne$X6),
                                  strand="*",
                                  name=paste0(cne$X1, ":", 
                                              (cne$X2+1), "-", cne$X3),
                                  itemRgb=chr2colour(cne$X1))
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
### makeBedWiggle: make the Ancora downloads-like bed files and wiggle files
###   from a GrangePairs of CNEs.
### 
makeBedWiggle <- function(x, outputDir=".",
                          firstGenome="first", secondGenome="second",
                          threshold="50_50"){
  if(!is(x, "GRangePairs")){
    stop("`x` must be a GRangePairs object!")
  }
  if(seqlengthsNA(x)){
    stop("seqlengths must be provided in `x`!")
  }
  
  ## make the bed files
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
                        name=paste(firstGenome, "CNEs", threshold),
                        description=paste(firstGenome, "CNEs", threshold)
                        )
  secondTrackLine <- new("BasicTrackLine", itemRgb=TRUE,
                         name=paste(secondGenome, "CNEs", threshold),
                         description=paste(secondGenome, "CNEs", threshold)
                         )
  export.bed(bedFirst, con=file.path(outputDir, 
                                     paste0("CNE_", firstGenome, "_",
                                            secondGenome, "_",
                                            threshold, ".bed")),
             trackLine=firstTrackLine)
  export.bed(bedSecond, con=file.path(outputDir, 
                                      paste0("CNE_", secondGenome, "_",
                                             firstGenome, "_",
                                             threshold, ".bed")),
             trackLine=secondTrackLine)
  
  # Make the wiggle files
  
}