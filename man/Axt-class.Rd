\name{Axt-class}
\Rdversion{1.1}
\docType{class}
\alias{Axt}
\alias{Axt-class}
\alias{[,Axt,ANY,ANY-method}
\alias{c,Axt-method}
\alias{length,Axt-method}
\alias{queryRanges}
\alias{queryRanges,Axt-method}
\alias{querySeqs}
\alias{querySeqs,Axt-method}
\alias{score,Axt-method}
\alias{symCount}
\alias{symCount,Axt-method}
\alias{targetRanges}
\alias{targetRanges,Axt-method}
\alias{targetSeqs}
\alias{targetSeqs,Axt-method}

\title{Class \code{"Axt"}}
\description{
  The Axt S4 object to hold a axt file.
}

\usage{
## Constructors:
Axt(targetRanges=GRanges(), targetSeqs=DNAStringSet(),
    queryRanges=GRanges(), querySeqs=DNAStringSet(),
    score=integer(0), symCount=integer(0), names=NULL)

## Accessor-like methods:
\S4method{targetRanges}{Axt}(x)
\S4method{targetSeqs}{Axt}(x)
\S4method{queryRanges}{Axt}(x)
\S4method{querySeqs}{Axt}(x)
\S4method{score}{Axt}(x)
\S4method{symCount}{Axt}(x)
## ... and more (see Methods)
}


\arguments{
  \item{targetRanges}{Object of class \code{"GRanges"}: 
    The ranges of net alignments on reference genome.}
  \item{targetSeqs}{Object of class \code{"DNAStringSet"}: 
    The alignment sequences of reference genome.}
  \item{queryRanges}{Object of class \code{"GRanges"}: 
    The ranges of net alignments on query genome.}
  \item{querySeqs}{Object of class \code{"DNAStringSet"}: 
    The alignment sequences of query genome.}
  \item{score}{Object of class \code{"integer"}: 
      The alignment score.}
  \item{symCount}{Object of class \code{"integer"}: 
    The alignment length.}
    \item{names}{\code{character}(): the names of axt alignments.}
  \item{x}{Object of class \code{"Axt"}:
    A Axt object.}
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "Axt", i = "ANY", j = "ANY")}: Axt getter}
    \item{c}{\code{signature(x = "Axt")}: Axt concatenator.}
    \item{length}{\code{signature(x = "Axt")}: Get the number of alignments.}
    \item{queryRanges}{\code{signature(x = "Axt")}: 
      Get the ranges of query genome.}
    \item{querySeqs}{\code{signature(x = "Axt")}: 
      Get the alignment sequences of query genome.}
    \item{score}{\code{signature(x = "Axt")}: Get the alignment score.}
    \item{symCount}{\code{signature(x = "Axt")}: Get the alignment lengths.}
    \item{targetRanges}{\code{signature(x = "Axt")}: 
      Get the ranges of reference genome.}
    \item{targetSeqs}{\code{signature(x = "Axt")}: 
      Get the alignment sequences of reference genome.}
	 }
}
\details{
  In \sQuote{axt} files and \code{Axt} object, the \sQuote{targetRanges}
  also have the alignments on positive strands.
  However, the \sQuote{queryRanges} can have alignments on negative strands,
  and the coordinates are based on negative strands, which is quite
  different from the convention in Bioconductor.
  To convert the coordinates of alignments on the negative strand to  
  the positive strand, use \code{normaliseStrand}.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{readAxt}}
  \code{\link{writeAxt}}
  \code{\link{subAxt}}
  \code{\link{fixCoordinates}}
  \code{\link{makeAxtTracks}}
}

\examples{
  library(GenomicRanges)
  library(Biostrings)
  ## Constructor
  targetRanges <- GRanges(seqnames=c("chr1", "chr1", "chr2", "chr3"),
                          ranges=IRanges(start=c(1, 20, 2, 3),
                                         end=c(10, 25, 10, 10)),
                          strand="+")
  targetSeqs <- DNAStringSet(c("ATTTTATGTG", "GGGAAG", "GGGCTTTTG",
                               "TTGTGTAG"))
  queryRanges <- GRanges(seqnames=c("chr1", "chr10", "chr10", "chr20"),
                         ranges=IRanges(start=c(1, 25, 50, 5),
                                        end=c(10, 30, 58, 12)),
                         strand="+")
  querySeqs <- DNAStringSet(c("ATTTAAAGTG", "GGAAAA", "GGGCTCTGG",
                              "TTAAATAA"))
  score <- c(246L, 4422L, 5679L, 1743L)
  symCount <- c(10L, 6L, 9L, 8L)
  axt <- Axt(targetRanges=targetRanges, targetSeqs=targetSeqs, 
             queryRanges=queryRanges, querySeqs=querySeqs, 
             score=score, symCount=symCount)
  
  ## getters
  names(axt)
  length(axt)
  first(axt)
  last(axt)
  seqnames(axt)
  strand(axt)
  seqinfo(axt)
  
  ## Vector methods
  axt[1]
  
  ## List methods
  unlist(axt)
  
  ## Combining
  c(axt, axt)
}

\keyword{classes}
