\name{cneMerge-methods}
\docType{methods}
\alias{cneMerge}
\alias{cneMerge,GRangePairs,GRangePairs-method}
\alias{cneMerge,CNE,missing-method}

\title{
  CNE merge function
}
\description{
  Removes the CNEs which overlap on both genomes.
}
\usage{
  cneMerge(cne12, cne21)
}

\arguments{
  \item{cne12}{
    A object of \code{CNE} or \code{GRangePairs}.
  }
  \item{cne21}{
    A object of \code{GRangePairs} object.
    When cne12 is a \code{CNE} object, cne21 can be missing.
  }
}


\value{
  A \code{GRangePairs} of CNEs or a \code{CNE} object is returned.
  In this table, the order of columns is consistent with cne1.
  For instance, if cne1 has the first three columns for zebrafish 
  and next three columns for human,
  in the merged table, the first three columns are 
  still the coordinates for zebrafish
  while the next three columns are the coordinates for human.
}

\author{
  Ge Tan
}




\examples{
  library(GenomicRanges)
  firstGRange <- GRanges(seqnames=c("chr1", "chr1", "chr2", "chr2", "chr5"),
                         ranges=IRanges(start=c(1, 20, 2, 3, 1),
                                        end=c(10, 25, 10, 10, 10)),
                         strand="+")
  lastGRange <- GRanges(seqnames=c("chr15", "chr10", "chr10", "chr10", "chr15"),
                        ranges=IRanges(start=c(1, 25, 50, 51, 5),
                                       end=c(8, 40, 55, 60, 10)),
                        strand="+")
  cne12 <- GRangePairs(firstGRange[1:3], lastGRange[1:3])
  cne21 <- GRangePairs(lastGRange[4:5], firstGRange[4:5])
  
  ## GRangePairs, GRangePairs
  cneMerge(cne12, cne21)
  
  ## CNE, missing
  cne <- CNE(assembly1Fn=file.path(system.file("extdata",
                                    package="BSgenome.Drerio.UCSC.danRer10"),
                        "single_sequences.2bit"),
             assembly2Fn=file.path(system.file("extdata",
                                    package="BSgenome.Hsapiens.UCSC.hg38"),
                        "single_sequences.2bit"), 
             window=50L, identity=50L,
             CNE12=cne12, CNE21=cne21, aligner="blat")
  cneMerge(cne)
}