## This is used to download the annotation from UCSC and 
### incooperate it for the tracks display

#.SUPPORTED_UCSC_TABLES = c(
#  ## tablename (unique key)   track             subtrack    auxiliary tablename
#  "knownGene",              "UCSC Genes",       NA,         
#  "refGene",                "RefSeq Genes",     NA,
#  "ensGene",                "Ensembl Genes",    NA
#  )

.SUPPORTED_UCSC_TABLES = list(
    "knownGene"    = c("knownGene", "kgXref"),
    "refGene"  = c("refGene"),
    "ensGene" = c("ensGene", "ensemblToGeneName")
    )

supportedUCSCtables = function(){
  .SUPPORTED_UCSC_TABLES
}

queryrefGene = function(con){
  query = "SELECT distinct name, name2, chrom, strand, exonStarts, exonEnds
                   FROM refGene
                   ORDER BY name2, name"
  ans = dbGetQuery(con, query)
  # process the ans into one exon per line
  exonStarts = strsplit(ans$exonStarts, ",")
  exonEnds = strsplit(ans$exonEnds, ",")
  stopifnot(all(sapply(exonStarts, length) == sapply(exonEnds, length)))
  repNum = sapply(exonStarts, length)
  res = data.frame(chromosome=rep(ans$chrom, repNum),
                   start=as.integer(unlist(exonStarts))+1,
                   end=as.integer(unlist(exonEnds)),
                   strand=rep(ans$strand, repNum),
                   gene=rep(ans$name2, repNum),
                   transcript=rep(ans$name, repNum),
                   symbol=rep(ans$name2, repNum))
  return(res)
}

queryknownGene = function(con){
  query = "SELECT distinct kgID, geneSymbol, chrom, strand, exonStarts, exonEnds
                   FROM knownGene, kgXref WHERE knownGene.name=kgXref.kgID 
                   ORDER BY geneSymbol, kgID"
  ans = dbGetQuery(con, query)
  # process the ans into one exon per line
  exonStarts = strsplit(ans$exonStarts, ",")
  exonEnds = strsplit(ans$exonEnds, ",")
  stopifnot(all(sapply(exonStarts, length) == sapply(exonEnds, length)))
  repNum = sapply(exonStarts, length)
  res = data.frame(chromosome=rep(ans$chrom, repNum),
                   start=as.integer(unlist(exonStarts))+1,
                   # The internal ucsc database use the 0-based start, 
                   # 1-based end. We only use 1-based.
                   end=as.integer(unlist(exonEnds)),
                   strand=rep(ans$strand, repNum),
                   gene=rep(ans$geneSymbol, repNum),
                   transcript=rep(ans$kgID, repNum),
                   symbol=rep(ans$geneSymbol, repNum))
  return(res)
}

queryensGene = function(con){
  query = "SELECT distinct chrom, strand, exonStarts, exonEnds, 
    ensGene.name2, ensGene.name, ensemblToGeneName.value
    FROM ensGene, ensemblToGeneName WHERE ensGene.name=ensemblToGeneName.name
    ORDER BY ensGene.name, ensemblToGeneName.value"
  ans = dbGetQuery(con, query)
  # process the ans into one exon per line
  exonStarts = strsplit(ans$exonStarts, ",")
  exonEnds = strsplit(ans$exonEnds, ",")
  stopifnot(all(sapply(exonStarts, length) == sapply(exonEnds, length)))
  repNum = sapply(exonStarts, length)
  res = data.frame(chromosome=rep(ans$chrom, repNum),
                   start=as.integer(unlist(exonStarts))+1,
                   end=as.integer(unlist(exonEnds)),
                   strand=rep(ans$strand, repNum),
                   gene=rep(ans$name2, repNum),
                   transcript=rep(ans$name, repNum),
                   symbol=rep(ans$value, repNum))
  return(res)
}

makeGeneDbFromUCSC = function(genome="hg19",
                              tablename="refGene",
                              host="genome-mysql.cse.ucsc.edu",
                              user="genome",
                              password=NULL,
                              dbnameSQLite="geneAnnotation.sqlite",
                              tablenameSQLite=paste(genome, tablename, sep="_"),
                              overwrite=FALSE 
                              ){
  if(!isSingleString(genome))
    stop("'genome' must be a single string")
  if(!isSingleString(tablename))
    stop("'tablename' must be a single string")
  if(!tablename %in% names(.SUPPORTED_UCSC_TABLES))
    stop("table \"", tablename, "\" is not supported")
  if(!isSingleString(host))
    stop("'url' must be a single string")
  con = dbConnect(MySQL(), user=user, password=password, 
                  dbname=genome, host=host)
  tableNames = .SUPPORTED_UCSC_TABLES[[tablename]] 
  message("Download the ", tablename, " table ... ")
  ans = switch(tablename,
               "refGene"=queryrefGene(con),
               "knownGene"=queryknownGene(con),
               "ensGene"=queryensGene(con)
               )
  dbDisconnect(con)
  # add the bin column
  ans$bin = binFromCoordRange(ans$start, ans$end)
  # reorder the columns, not necessary
  ans = ans[ ,c("bin","chromosome","start","end","strand",
                "gene", "transcript","symbol")]
  con = dbConnect(SQLite(), dbname=dbnameSQLite)
  dbWriteTable(con, tablenameSQLite, ans, overwrite=overwrite)
  dbDisconnect(con)
}



