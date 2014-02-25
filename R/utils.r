### -----------------------------------------------------------------
### Compare two DNAStringSet. Match is TRUE, Mismatch is FALSE.
### Not used, exported so far.
compDNAStringSet <- function(DNAStringSet1, DNAStringSet2){
  tmp <- cbind(strsplit(as.character(DNAStringSet1), ""), 
               strsplit(as.character(DNAStringSet2), ""))
  apply(tmp, 1, function(x){x[[1]] == x[[2]]})
}
#system.time(foo<-compDNAStringSet(targetSeqs(myAxt), querySeqs(myAxt)))
#system.time(foo1<-RleList(foo))


### -----------------------------------------------------------------
### For a GRanges object of filter, make the revered GRanges. 
### chromSize needs to be known.
### Not used, exported so far.
makeReversedFilter <- function(qFilter, chromSizes){
  revFilterBed <- GRanges(seqnames=seqnames(qFilter),
                          ranges=IRanges(start=chromSizes[
                                         as.character(seqnames(qFilter))] - 
                                         end(qFilter),
                                         end=chromSizes[
                                         as.character(seqnames(qFilter))] - 
                                         start(qFilter)
                                         ),
                          strand=Rle("-"))
  return(revFilterBed)
}

### -----------------------------------------------------------------
### Generate a translation from sequence index to alignment index.
### Not used, exported so far.
seqToAlignment <- function(DNAStringSet){
  foo <- strsplit(as.character(DNAStringSet), "")
  foo <- lapply(foo, function(x){grep("-", x, invert=TRUE)})
  return(foo)
}

### -----------------------------------------------------------------
### rever the cigar string. i.e. 20M15I10D will be reversed to 10D15I20M.
### EXPORTED!
reverseCigar <- function(cigar, ops=CIGAR_OPS){
  #cigar = sapply(splitCigar(cigar), function(x){
  #               paste0(rev(x[[2]]), rev(rawToChar(x[[1]], multiple=TRUE)), 
  #                      collapse="")
  #                       }
  # splitCigar is deprecated...Before I am in the Bioconductor..
  # some new cigar utilities functions.
  #)
  cigarOps <- lapply(explodeCigarOps(cigar, ops=ops), rev)
  cigarOpsLengths <- lapply(explodeCigarOpLengths(cigar, ops=ops), rev)
  cigar <- mapply(paste0, cigarOpsLengths, cigarOps, collapse="")
  return(cigar)
}

### -----------------------------------------------------------------
### Better system interface
### Not exported.
my.system <- function(cmd, echo=TRUE, intern=FALSE, ...){
  if (echo){
    message(cmd)
  }
  res <- system(cmd, intern=intern, ...)
  if (!intern){
    stopifnot(res == 0)
  }
  return(res)
}


### -----------------------------------------------------------------
### Return the bin number that should be assigned to 
### a feature spanning the given range. * USE THIS WHEN CREATING A DB *
### Exported!
.validateBinRanges <- function(starts, ends){
  if(any(ends <= 0 | starts <= 0)){
    stop("starts and ends must be positive integers!")
  }
  if(any(ends >= 2^30 | starts >= 2^30)){
    stop("starts and ends out of range (max is 2Gb)")
  }
  if(any(starts > ends)){
    stop("starts must be equal or smaller than ends!")
  }
  return(TRUE) 
}

binFromCoordRange <- function(starts, ends){
  .validateBinRanges(starts, ends)
  bins <- .Call2("bin_from_coord_range", as.integer(starts), 
                 as.integer(ends), PACKAGE="CNEr")
  return(bins)
}

### -----------------------------------------------------------------
### Return the set of bin ranges that overlap a given coordinate range. 
### It is usually more convenient to use bin_restriction string 
### than to use this method directly.
### EXPORTED!
binRangesFromCoordRange <- function(start, end){
  stopifnot(length(start) == 1 && length(end) == 1)
  .validateBinRanges(start, end)
  binRanges <- .Call2("bin_ranges_from_coord_range", 
                      as.integer(start), as.integer(end), PACKAGE="CNEr")
  return(binRanges)
}

### -----------------------------------------------------------------
### Given a coordinate range, return a string to be used in the WHERE 
### section of a SQL SELECT statement that is to 
### select features overlapping a certain range. * USE THIS WHEN QUERYING A DB *
### EXPORTED!
binRestrictionString <- function(start, end, field="bin"){
  binRanges <- binRangesFromCoordRange(start, end)
  cmdString <- mapply(function(x,y, field){
                      if(x==y){
                        paste(field, "=", x)
                      }else{
                        paste(field, ">=", x, "and", field, "<=", y)
                      }
                     }, binRanges[ ,1], binRanges[ ,2], field=field
                     )
  cmdString <- paste(cmdString, collapse=") or (")
  cmdString <- paste0("((", cmdString, "))")
  return(cmdString)
}





#get_cne_ranges_in_region = function(CNE, whichAssembly=c(1,2), 
#                                    chr, CNEstart, CNEend, min_length){
#  ## This CNE data.frame does not have the bin column yet. 
#  ## I am not sure whether it is necessary to add this column in R since 
#  ## it's quiet fast to select the cnes which meet the criteria (~0.005 second).
#  if(whichAssembly == 1)
#    res = subset(CNE, chr1==chr & start1>=CNEstart & end1<=CNEend & 
#                 end1-start1+1>=min_length & end2-start2+1>=min_length, 
#                 select=c("start1", "end1"))
#  else if(whichAssembly == 2)
#    res = subset(CNE, chr2==chr & start2>=CNEstart & end2<=CNEend & 
#                 end1-start1+1>=min_length & end2-start2+1>=min_length, 
#                 select=c("start2", "end2"))
#  else
#    stop("whichAssembly should be 1 or 2")
#  # Here we return a IRanges object to store the start and end
#  res = IRanges(start=res[ ,1], end=res[ ,2])
#  return(res)
#}
#
#
#get_cne_ranges_in_region_partitioned_by_other_chr = 
#  function(CNE, whichAssembly=c(1,2), chr, CNEstart, CNEend, min_length){
#  if(whichAssembly == 1)
#    res = subset(CNE, chr1==chr & start1>=CNEstart & end1<=CNEend & 
#                 end1-start1+1>=min_length & end2-start2+1>=min_length, 
#                 select=c("chr2", "start1", "end1"))
#  else if(whichAssembly == 1)
#    res = subset(CNE, chr2==chr & start2>=CNEstart & end2<=CNEend & 
#                 end1-start1+1>=min_length & end2-start2+1>=min_length, 
#                 select=c("chr1", "start2", "end2"))
#  else
#    stop("whichAssembly should be 1 or 2")
#  # Here we return a GRanges object.
#  res = GRanges(seqnames=res[ ,1], 
#                ranges=IRanges(start=res[ ,2], end=res[ ,3]))
#  return(res)
#}

### -----------------------------------------------------------------
### save the CNE tables into a local SQLite database
### Exported!!
setMethod("saveCNEToSQLite",
          signature(CNE="data.frame", tableName="character"),
          function(CNE, dbName, tableName, overwrite=FALSE){
            ## tableName should be in the format "danRer7_hg19_49_50"
            if(!grepl("^.+_.+_\\d+_\\d+$", tableName))
              stop("The tableName should be in the format danRer7_hg19_49_50.")
            CNE$bin1 <- binFromCoordRange(CNE$start1, CNE$end1)
            CNE$bin2 <- binFromCoordRange(CNE$start2, CNE$end2)
            # reorder it
            CNE <- CNE[ ,c("bin1", "chr1", "start1", "end1", "bin2",
                           "chr2", "start2", "end2", "strand", 
                           "similarity", "cigar")]
            con <- dbConnect(SQLite(), dbname=dbName)
            on.exit(dbDisconnect(con))
            dbWriteTable(con, tableName, CNE, row.names=FALSE, overwrite=overwrite)
          }
          )

setMethod("saveCNEToSQLite",
          signature(CNE="CNE", tableName="missing"),
          function(CNE, dbName, tableName, overwrite=FALSE){
            tableNames = names(CNE@CNERepeatsFiltered)
            for(i in 1:length(CNE@CNERepeatsFiltered)){
              saveCNEToSQLite(CNE@CNERepeatsFiltered[[i]],
                              dbName=dbName, tableName=tableNames[i],
                              overwrite=overwrite)
            }
          }
          )
          



### -----------------------------------------------------------------
### read CNE ranges from a local SQLite database.
### Exported!
readCNERangesFromSQLite <- function(dbName, tableName, chr, start, end, 
                                    whichAssembly=c("L","R"), minLength=NULL){
  nrGraphs <- 1
  ## Let's make nrGraphs=1, make all the cnes together.
  if(!is(start, "integer"))
    stop("start must be an integer!")
  if(!is(end, "integer"))
    stop("end must be an integer!")
  CNEstart <- start
  CNEend <- end
  whichAssembly <- match.arg(whichAssembly)
  con <- dbConnect(SQLite(), dbname=dbName)
  on.exit(dbDisconnect(con))
  if(nrGraphs == 1){
    sqlCmd <- switch(whichAssembly,
                     "L"=paste("SELECT start1,end1 from", tableName, 
                               "WHERE chr1=", paste0("'", chr, "'"), 
                               "AND start1 >=", CNEstart, "AND end1 <=", 
                               CNEend, "AND", 
                               binRestrictionString(CNEstart, CNEend, "bin1")),
                     "R"=paste("SELECT start2,end2 from", tableName, 
                               "WHERE chr2=", paste0("'", chr, "'"), 
                               "AND start2 >=", CNEstart, "AND end2 <=", 
                               CNEend, "AND", 
                               binRestrictionString(CNEstart, CNEend, "bin2"))
                     )
    if(!is.null(minLength))
      sqlCmd <- paste(sqlCmd, "AND end1-start1+1 >=", minLength, 
                      "AND end2-start2+1 >=", minLength)
    fetchedCNE <- dbGetQuery(con, sqlCmd)
    fetchedCNE <- IRanges(start=fetchedCNE[ ,1], end=fetchedCNE[, 2])
  }else if(nrGraphs > 1){
    sqlCmd <- switch(whichAssembly,
                     "L"=paste("SELECT chr2,start1,end1 from", tableName, 
                               "WHERE chr1=", paste0("'", chr, "'"), 
                               "AND start1 >=", CNEstart, "AND end1 <=", 
                               CNEend, "AND", 
                               binRestrictionString(CNEstart, CNEend, "bin1")),
                     "R"=paste("SELECT chr1,start2,end2 from", tableName, 
                               "WHERE chr2=", paste0("'", chr, "'"), 
                               "AND start2 >=", CNEstart, "AND end2 <=", 
                               CNEend, "AND", 
                               binRestrictionString(CNEstart, CNEend, "bin2"))
                     )
    if(!is.null(minLength))
      sqlCmd <- paste(sqlCmd, "AND end1-start1+1 >=", minLength, 
                      "AND end2-start2+1 >=", minLength)
    fetchedCNE <- dbGetQuery(con, sqlCmd)
    fetchedCNE <- GRanges(seqnames=fetchedCNE[ ,1], 
                          ranges=IRanges(start=fetchedCNE[ ,2], 
                                         end=fetchedCNE[ ,3]))
  }
  return(fetchedCNE)
}

queryAnnotationSQLite <- function(dbname, tablename, chr, start, end){
  con <- dbConnect(SQLite(), dbname=dbname)
  query <- paste("SELECT * from", tablename, "WHERE", 
                binRestrictionString(start, end, "bin"), "AND", 
                "chromosome=", paste0("'", chr, "'"), 
                "AND start >=", start, "AND end <=", end)
  ans <- dbGetQuery(con, query)
  ans <- ans[ ,c("chromosome", "start", "end", "strand", 
                "transcript", "gene", "symbol")]
}

###------------------------------------------------------------------
### fetchChromSizes fetches the chromosome sizes.
### Exported!
fetchChromSizes <- function(assembly){
  # UCSC
  message("Trying UCSC...")
  con <- try(dbConnect(MySQL(), user="genome", password="", 
                       dbname=assembly, host="genome-mysql.cse.ucsc.edu"), 
             silent=TRUE)
  if(class(con) != "try-error"){
    on.exit(dbDisconnect(con))
    sqlCmd <- "SELECT chrom,size FROM chromInfo ORDER BY size DESC"
    ans <- try(dbGetQuery(con, sqlCmd))
    if(class(ans) == "try-error"){
      return(NULL)
    }else{
      ans <- Seqinfo(seqnames=ans$chrom, seqlengths=ans$size, genome=assembly)
      return(ans)
    }
  }
  # other sources? Add later
  return(NULL)
}


