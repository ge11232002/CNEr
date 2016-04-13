\name{subAxt-methods}
\docType{methods}
\alias{subAxt}
\alias{subAxt-methods}
\alias{subAxt,Axt,character,missing,missing-method}
\alias{subAxt,Axt,character,numeric,numeric-method}
\alias{subAxt,Axt,character,integer,integer-method}
\alias{subAxt,Axt,GRanges,missing,missing-method}

\title{Subset an \code{Axt} object}
\description{
  A 'subAxt' method for extracting a set of alignments from
  an \code{Axt} object.
}

\usage{
  subAxt(x, chr, start, end, select=c("target", "query"), qSize=NULL)
}

\arguments{
  \item{x}{
    A object of \code{Axt}.
  }
  \item{chr}{
    A object of \code{character} containing the names of the sequences 
    in 'x' where to get the alignments from, or a \code{GRanges} object
    when 'start' and 'end' are missing.
  }
  \item{start, end}{
    A object of \code{integer}() or missing.
    These ranges should be based on the positive strand.
    When select is "query",
    the reverse complement alignments which lay inside this range will also be 
    selected.
  }
  \item{select}{
    When select is "target", 
    the subset criteria is for target alignments in \code{Axt}.
    When select is "query", 
    the subset criteria is for query alignments in \code{Axt}.
  }
  \item{qSize}{
    When select is "query", qSize must be provided as the same length of 'chr'.
    It holds the chromosome length of 'chr'.
  }
}

\details{
  Usually when we want to subset some axts from a \code{Axt} object,
  we care about all the axts within certain range.
  The axts can come from the axt file with chr as reference 
  (i.e., target sequence),
  or the axt file with chr as query sequence.
  When the chr is query sequence, it can be on the negative strand.
  Hence, the size of chromosome is necessary to 
  convert the search range to a range on negative strand coordinate.

  When one Axt alignment is partially overlapped with the range, 
  the full length of the Axt alignment will be extracted.
}

\value{
  An extracted \code{Axt} object is returned.
}

\author{
  Ge Tan
}

\examples{
  library(GenomicRanges)
  ## Prepare the axt object
  axtFilesHg19DanRer7 <- file.path(system.file("extdata", package="CNEr"), 
                                  "hg19.danRer7.net.axt")
  axtHg19DanRer7 <- readAxt(axtFilesHg19DanRer7)
  
  ## "character", "integer", "integer" on "target" sequence
  subAxt(axtHg19DanRer7, chr="chr11", start=31500000L, end=32500000L,
         select="target")
  
  ## "GRanges"" on "target" sequence
  searchGRanges <- GRanges(seqnames="chr11",
                           ranges=IRanges(start=31500000,
                                          end=32500000),
                           strand="+")
  subAxt(axtHg19DanRer7, searchGRanges, select="target")
  
  ## multiple "character", "integer", "integer" on "target" sequence
  subAxt(axtHg19DanRer7, chr="chr11", start=c(31082021L, 32461267L),
         end=c(31082862L, 32461581L), select="target")
  
  ## "character" only on "target" sequence
  subAxt(axtHg19DanRer7, chr="chr11", select="target")
  
  ## GRanges on "query" sequence
  searchGRanges <- GRanges(seqnames="chr25",
                           ranges=IRanges(start=15559655,
                                          end=15575192),
                           strand="+")
  subAxt(axtHg19DanRer7, searchGRanges, select="query", 
         qSize=c("chr25"=38499472L))
}