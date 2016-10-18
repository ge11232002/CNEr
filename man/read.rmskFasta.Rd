\name{read.rmskFasta}
\alias{read.rmskFasta}
\title{
  Read a soft repeat masked fasta
}
\description{
  Read a soft repeat masked fasta file into a \code{GRanges} object.
}
\usage{
  read.rmskFasta(fn)
}

\arguments{
  \item{fn}{
    \code{character}(1): The filename of the soft repeat masked fasta.
  }
}
\details{
  Only the lower case based ("a", "c", "g", "t") are considered in the
  soft repeat masked fasta.
}
\value{
  \code{GRanges} object with coordinates of repeat masked regions.
}
\author{
  Ge Tan
}
\seealso{
  \code{\link{read.rmMask.GRanges}}
}
\examples{
  fn <- file.path(system.file("extdata", package="CNEr"),
                  "rmsk.fa")
  read.rmskFasta(fn)
}