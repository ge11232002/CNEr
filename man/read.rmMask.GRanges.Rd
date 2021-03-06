\name{read.rmMask.GRanges}
\alias{read.rmMask.GRanges}

\title{
  Read a RepeatMasker .out file 
}
\description{
  Read a RepeatMasker .out file into a \code{GRanges} object.
}
\usage{
  read.rmMask.GRanges(fn)
}

\arguments{
  \item{fn}{
    \code{character}(1): the filename of a RepeatMasker .out file.
  }
}

\value{
  A \code{GRanges} object with metadata columns containing the
  name of the matching interspersed repeat,
  the class of the repeat and the Smith-Waterman score of the match.
}
\references{
  \url{http://www.repeatmasker.org/webrepeatmaskerhelp.html}
}
\author{
  Ge Tan
}

\examples{
  fn <- system.file("extdata", "ce2chrM.fa.out", package="IRanges")
  read.rmMask.GRanges(fn)
}
