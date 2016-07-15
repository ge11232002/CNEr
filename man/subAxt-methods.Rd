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
  A \sQuote{subAxt} method for extracting a set of alignments from
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
    In the case of \code{GRanges}, the strand information is ignored.
  }
  \item{start, end}{
    A object of \code{integer}() or missing.
    These ranges should be based on the positive strand.
    When select is "query",
    the reverse complement alignments which lay inside this range will also be 
    selected.
  }
  \item{select}{
    When select is \sQuote{target},
    the subset criteria is for target alignments in \code{Axt}.
    When select is \sQuote{query}, 
    the subset criteria is for query alignments in \code{Axt}.
  }
  \item{qSize}{
    \code{integer}(n): When select is \sQuote{query}, qSize must exist in 
      \sQuote{x} or is provided as vector of chromosome lengths.
  }
}

\details{
  Usually when we want to subset some axts from a \code{Axt} object,
  we care about all the axts within certain range.
  The axts can come from the axt file with chr as reference 
  (\emph{i.e.}, target sequence),
  or the axt file with chr as query sequence.
  When the chr is query sequence, it can be on the negative strand.
  Hence, the size of chromosome is necessary to 
  convert the search range to a range on negative strand coordinate.

  When one Axt alignment is partially overlapped with the range, 
  the whole Axt alignment will be extracted.
}

\value{
  An extracted \code{Axt} object is returned.
}

\author{
  Ge Tan
}

\examples{
  library(GenomicRanges)
  library(rtracklayer)
  ## Prepare the axt object
  axtFilesHg38DanRer10 <- file.path(system.file("extdata", package="CNEr"), 
                                  "hg38.danRer10.net.axt")
  axtHg38DanRer10 <- readAxt(axtFilesHg38DanRer10)
  
  ## "character", "integer", "integer" on "target" sequence
  subAxt(axtHg38DanRer10, chr="chr1", start=148165963L, end=222131835L,
         select="target")
  
  ## "GRanges"" on "target" sequence
  searchGRanges <- GRanges(seqnames="chr1",
                           ranges=IRanges(start=148165963L,
                                          end=222131835L),
                           strand="+")
  subAxt(axtHg38DanRer10, searchGRanges, select="target")
  
  ## multiple "character", "integer", "integer" on "target" sequence
  subAxt(axtHg38DanRer10, chr=c("chr1", "chr13"), 
         start=c(148165963L, 94750629L),
         end=c(222131835L, 94966991L), select="target")
  
  ## "character" only on "target" sequence
  subAxt(axtHg38DanRer10, chr="chr1", select="target")
  
  ## GRanges on "query" sequence
  danRer10Seqlengths <- seqlengths(TwoBitFile(file.path(system.file("extdata",
                                   package="BSgenome.Drerio.UCSC.danRer10"),
                                   "single_sequences.2bit")))
  searchGRanges <- GRanges(seqnames="chr6",
                           ranges=IRanges(start=25825774,
                                          end=26745499),
                           strand="+")
  subAxt(axtHg38DanRer10, searchGRanges, select="query", 
         qSize=danRer10Seqlengths["chr6"])
}