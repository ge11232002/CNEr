\name{CNEDensity-methods}
\docType{methods}
\alias{CNEDensity}
\alias{CNEDensity-methods}
\alias{CNEDensity,ANY,character,character,missing,missing-method}
\alias{CNEDensity,ANY,missing,character,character,character-method}

\title{
  CNEDensity function
}

\description{
  This function queries the database and generates the CNEs' density values.
}

\usage{
CNEDensity(dbName, tableName, chr, start, end, 
           whichAssembly=c("first", "second"),
           windowSize=300, minLength=NULL)
}

\arguments{
  \item{dbName}{
    \code{character}(1): the path of the local SQLite database.
  }
  \item{tableName}{
    \code{character}(1): the name of table for this CNE data table.
    It can be missing when assembly1, assembly2 and threshold are provided.
  }
  \item{chr}{
     \code{character}(1): the chromosome to query.
  }
  \item{start, end}{
    \code{integer}(1): the start and end coordinate to fetch the CNEs.
  }
  \item{whichAssembly}{
    \code{character}(1):
    the genome to fetch is in the \sQuote{first} column or \sQuote{second}
    column of the table.
  }
  \item{windowSize}{
    \code{integer}(1): the window size in kb that is used to smooth the CNEs.
  }
  \item{minLength}{
    \code{integer}(1): the minimal length of CNEs to fetch.
  }
}

\section{Methods}{
\describe{

\item{\code{signature(tableName = "character",
                      assembly1 = "character", assembly2 = "missing",
                      threshold = "missing")}}{
}

\item{\code{signature(tableName = "missing", assembly1 = "character",
                      assembly2 = "character", threshold = "character")}}{
}
}}

\value{
  A \code{GRanges} object with density values is returned. 
}

\author{
  Ge Tan
}


\examples{
dbName <- file.path(system.file("extdata", package="CNEr"),
                    "danRer10CNE.sqlite")
genome <- "danRer10"
chr <- "chr6"
start <- 24000000L
end <-  27000000L
windowSize <- 200L
minLength <- 50L
cneDanRer10Hg38_45_50 <- 
  CNEDensity(dbName=dbName, 
             tableName="danRer10_hg38_45_50", 
             whichAssembly="first", chr=chr, start=start,
             end=end, windowSize=windowSize, 
             minLength=minLength)
cneDanRer10Hg38_49_50 <-
  CNEDensity(dbName=dbName, 
             tableName="danRer10_hg38_49_50", 
             whichAssembly="first", chr=chr, start=start,
             end=end, windowSize=windowSize, 
             minLength=minLength)
}