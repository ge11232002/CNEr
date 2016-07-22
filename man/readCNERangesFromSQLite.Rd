\name{readCNERangesFromSQLite}
\alias{readCNERangesFromSQLite}

\title{
  readCNERangesFromSQLite function
}
\description{
Query the SQLite database based on chromosome, 
coordiantes and some other criterias.
Usually not to be used directly.
For the CNE density plot, \code{fetchCNEDensity} function should be used.
}
\usage{
readCNERangesFromSQLite(dbName, tableName, chr=NULL, start=NULL, end=NULL,
                        whichAssembly=c("first","second"), minLength=NULL,
                        tAssemblyFn=NULL, qAssemblyFn=NULL)
}

\arguments{
  \item{dbName}{
    A object of \code{character}, the path of the local SQLite database.
  }
  \item{tableName}{
    A object of \code{character}, the name of table for this CNE data table.
  }
  \item{chr}{
    \code{character}(n): the chromosomes to query.
    When it's NULL, all CNEs will be returned.
  }
  \item{start, end}{
    \code{integer}(n): the start and end coordiantes to fetch the CNEs.
  }
  \item{whichAssembly}{
    \code{character}(1):
    The coordinates to fetch CNEs are based on \sQuote{first} genome or 
    \sQuote{last} genome.
  }
  \item{minLength}{
    \code{integer}(1): the minimal length for selected CNEs.
    The pair of CNEs both must be longer than this \code{minLength}.
  }
  \item{tAssemblyFn,qAssemblyFn}{
    \code{character}(1): filename of the \sQuote{twoBit} or \sQuote{fasta} 
    file for the target and query genome.
  }
}

\value{
  A object of \code{GRangePairs} is returned.
}

\author{
  Ge Tan
}

\examples{
  dbName <- file.path(system.file("extdata", package="CNEr"),
                    "danRer10CNE.sqlite")
  tableName <- "danRer10_hg38_45_50"
  
  qAssemblyFn <- file.path(system.file("extdata",
                           package="BSgenome.Hsapiens.UCSC.hg38"),
                           "single_sequences.2bit")
  tAssemblyFn <- file.path(system.file("extdata",
                           package="BSgenome.Drerio.UCSC.danRer10"),
                           "single_sequences.2bit")
                           
  ## single chr, start, end
  chr <- "chr6"
  start <- 24000000L
  end <-  27000000
  minLength <- 50L
  fetchedCNERanges <- readCNERangesFromSQLite(dbName, tableName, chr, 
                                              start, end, 
                                              whichAssembly="first",
                                              minLength=minLength,
                                              tAssemblyFn=tAssemblyFn,
                                              qAssemblyFn=qAssemblyFn)
  ## multiple chr, start, end
  chr=c("chr1", "chr3")
  start=c(90730248, 137523122)
  end=c(90730300, 137523190)
  fetchedCNERanges <- readCNERangesFromSQLite(dbName, tableName, chr, 
                                              start, end, 
                                              whichAssembly="second",
                                              minLength=minLength)
  
  ## chr, NULL, NULL
  fetchedCNERanges <- readCNERangesFromSQLite(dbName, tableName, chr, 
                                              start=NULL, end=NULL,
                                              whichAssembly="second",
                                              minLength=minLength)
}