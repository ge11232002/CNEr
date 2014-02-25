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
This function queries the database and generates the CNEs density values.
}

\usage{
CNEDensity(dbName, tableName, assembly1, assembly2, threshold,
           chr, start, end, windowSize, minLength=NULL)
}

\arguments{
  \item{dbName}{
    A object of \code{character}, the path of the local SQLite database.
  }
  \item{tableName}{
    A object of \code{character}, the name of table for this CNE data table.
    It can be missing when assembly1, assembly2 and threshold are provided.
  }
  \item{assembly1}{
    A object of \code{character}, the assembly to search.
  }
  \item{assembly2}{
    The comparison assembly. It can be missing when tableName is provided.
  } 
  \item{threshold}{
    The threshold to search. It can be missing when tableName is provided.
  }
  \item{chr}{
     A object of \code{character}, the chromosome to query.
  }
  \item{start, end}{
    A object of \code{integer}, the start and end coordiante to fetch the CNEs.
  }
  \item{windowSize}{
    A object of \code{integer}, the window size in kb used to smooth the CNEs.
  }
  \item{minLength}{
    A object of \code{integer}, the minimal length of CNEs to fetch.
  }
}

\section{Methods}{
\describe{

\item{\code{signature(tableName = "character", assembly1 = "character", assembly2 = "missing", threshold = "missing")}}{
}

\item{\code{signature(tableName = "missing", assembly1 = "character", assembly2 = "character", threshold = "character")}}{
}
}}

\value{
  A \code{matrix} is returned. 
  The first column is the coordinates and the second column 
  is the density values.
}

\author{
Ge Tan
}


\examples{
    dbName <- file.path(system.file("extdata", package="CNEr"),
                        "cne.sqlite")
    chr <- "chr11"
    start <- 31000000L
    end <-  33000000L
    windowSize <- 300L
    minLength <- 50L
    cneHg19DanRer7_45_50 <- 
      CNEDensity(dbName=dbName, 
                 tableName="danRer7_hg19_45_50", 
                 assembly1="hg19", chr=chr, start=start,
                 end=end, windowSize=windowSize, 
                 minLength=minLength)
    cneHg19DanRer7_48_50 <- 
      CNEDensity(dbName=dbName, 
                 tableName="danRer7_hg19_45_50", 
                 assembly1="hg19", chr=chr, start=start,
                 end=end, windowSize=windowSize, 
                 minLength=minLength)
    cneHg19DanRer7_49_50 <- 
      CNEDensity(dbName=dbName, 
                 tableName="danRer7_hg19_45_50", 
                 assembly1="hg19", chr=chr, start=start,
                 end=end, windowSize=windowSize, 
                 minLength=minLength)
}


