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
                        whichAssembly=c("first","last"), minLength=NULL)
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
    The coordinates to fetch CNEs are based on "first" genome or "last" genome.
  }
  \item{minLength}{
    \code{integer}(1): the minimal length for selected CNEs.
    The pair of CNEs both must be longer than this \code{minLength}.
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
                      "cne.sqlite")
  tableName <- "danRer10_hg38_45_50"
  
  ## single chr, start, end
  chr <- "chr6"
  start <- 24000000L
  end <-  27000000
  minLength <- 50L
  fetchedCNERanges <- readCNERangesFromSQLite(dbName, tableName, chr, 
                                              start, end, 
                                              whichAssembly="first",
                                              minLength=minLength)
  ## multiple chr, start, end
  chr=c("chr1", "chr3")
  start=c(90730248, 137523122)
  end=c(90730300, 137523190)
  fetchedCNERanges <- readCNERangesFromSQLite(dbName, tableName, chr, 
                                              start, end, 
                                              whichAssembly="last",
                                              minLength=minLength)
  
  ## chr, NULL, NULL
  fetchedCNERanges <- readCNERangesFromSQLite(dbName, tableName, chr, 
                                              start=NULL, end=NULL,
                                              whichAssembly="last",
                                              minLength=minLength)
}