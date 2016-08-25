\name{fixCoordinates}
\alias{fixCoordinates}
\alias{fixCoordinates,Axt-method}
\title{
  Fix the coordinates in \code{Axt} object
}
\description{
  In \sQuote{axt} file and \code{Axt} object, the coordinates of negative
  query alignments are relative to the reverse-complemented coordinates of 
  its chromosome.
  This is different from the convention in Bioconductor.
  This function fixes the coordinates which are always relative to the
  positive strand.
}
\usage{
  fixCoordinates(x)
}
\arguments{
  \item{x}{
    \code{Axt} object.
  }
}
\details{
  In \code{Axt}, the \sQuote{strand} is for the aligning organism.
  If the strand value is \dQuote{-}, the values of the aligning organism's 
  start and end fields are relative to the reverse-complemented coordinates of 
  its chromosome. 
}
\value{
  A \code{Axt} object.
}
\author{
  Ge Tan
}

\examples{
  axtFnDanRer10Hg38 <- file.path(system.file("extdata", package="CNEr"),
                                 "danRer10.hg38.net.axt")
  qAssemblyFn <- file.path(system.file("extdata",
                                       package="BSgenome.Hsapiens.UCSC.hg38"),
                           "single_sequences.2bit")
  tAssemblyFn <- file.path(system.file("extdata",
                                       package="BSgenome.Drerio.UCSC.danRer10"),
                           "single_sequences.2bit")
  axtDanRer10Hg38 <- readAxt(axtFnDanRer10Hg38, tAssemblyFn=tAssemblyFn,
                             qAssemblyFn=qAssemblyFn)
  ## Fix the coordinates
  fixCoordinates(axtDanRer10Hg38)
  
  ## Restore it
  fixCoordinates(fixCoordinates(axtDanRer10Hg38))
}