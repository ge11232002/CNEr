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
    score=integer(), symCount=integer())

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

\author{
  Ge Tan
}

\seealso{
  \code{\link{readAxt}}
  \code{\link{writeAxt}}
  \code{\link{subAxt}}
}

\examples{
showClass("Axt")
}

\keyword{classes}
