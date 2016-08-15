\name{psubAxt}
\alias{psubAxt}
\title{
  Parallel subset of Axt alignment
}
\description{
  Give two \code{GRanges} objects, select the \code{Axt} alignments
  whose the target and query alignments are both within each pair of ranges.
}
\usage{
  psubAxt(x, targetSearch, querySearch)
}
\arguments{
  \item{x}{
    \code{Axt} object.
  }
  \item{targetSearch,querySearch}{
    \code{GRanges} objects: the ranges to keep for target and query alignments.
    They must have same lengths and the strand information is ignored.
  }
}
\details{
  The \sQuote{targetSearch} and \sQuote{querySearch} have the coordinates
  relative to the positive strand.
  For each pair of the ranges, the alignments lie within both the target
  and query range are kept.
}
\value{
  A \code{Axt} object.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{psubAxt}}
}

\examples{
  library(GenomicRanges)
  tAssemblyFn <- file.path(system.file("extdata",
                             package="BSgenome.Drerio.UCSC.danRer10"),
                           "single_sequences.2bit")
  qAssemblyFn <- file.path(system.file("extdata",
                             package="BSgenome.Hsapiens.UCSC.hg38"),
                           "single_sequences.2bit")
  axtFn <- file.path(system.file("extdata", package="CNEr"), 
                     "danRer10.hg38.net.axt")
  axt <- readAxt(axtFn, tAssemblyFn, qAssemblyFn)

  targetSearch <- GRanges(seqnames=c("chr6"),
                          ranges=IRanges(start=c(24000000, 26900000),
                                         end=c(24060000, 26905000)),
                          strand="+"
                          )
  querySearch <- GRanges(seqnames=c("chr7", "chr2"),
                         ranges=IRanges(start=c(12577000, 241262700),
                                        end=c(12579000, 241268600)),
                         strand="+"
                         )
  psubAxt(axt, targetSearch, querySearch)
}