\name{makeAxtTracks}
\alias{makeAxtTracks}
\title{
  makeAxtTracks
}
\description{
  Make the bed tracks for the \sQuote{Axt} alignment.
}
\usage{
  makeAxtTracks(x)
}
\arguments{
  \item{x}{
    A \code{Axt} object.
  }
}
\details{
  The coordinates of query \sQuote{Axt} alignment are fixed to be relative to
  positive strand before output into \sQuote{bed} file.
}
\value{
  A list of GRanges for target and query alignments.
  The two output \sQuote{bed} files are 
  \dQuote{targetAxt.bed} and \dQuote{queryAxt.bed}.
}
\author{
  Ge Tan
}
\seealso{
  \code{\link{fixCoordinates}}
}
\examples{
  tAssemblyFn <- file.path(system.file("extdata",
                                       package="BSgenome.Drerio.UCSC.danRer10"),
                           "single_sequences.2bit")
  qAssemblyFn <- file.path(system.file("extdata",
                                       package="BSgenome.Hsapiens.UCSC.hg38"),
                           "single_sequences.2bit")
  axtFn <- file.path(system.file("extdata", package="CNEr"), 
                     "danRer10.hg38.net.axt")
  axt <- readAxt(axtFn, tAssemblyFn, qAssemblyFn)
  makeAxtTracks(axt)
}
