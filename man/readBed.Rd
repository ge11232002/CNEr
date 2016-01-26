\name{readBed}
\alias{readBed}
\title{
  readBed
}
\description{
  Read the coordiantes information from a bed file.
}
\usage{
  readBed(bedFile)
}
\arguments{
  \item{bedFile}{
    The \code{character}(1) file name of the 'bed' file to read.
  }
}
\details{
  This function is designed to read the bed file for the 
  first three columns,
  \emph{i.e.}, "chrom", "chromStart", "chromEnd".
  The strand information is also stored when available.
  
  In bed file, the "chromStart" is on 0-based coordinate 
  while "chromEnd" is on 1-based coordinate.
  For example, the first 100 bases of a chromosome are 
  defined as "chromStart"=0, "chromEnd"=100, 
  and span the bases numbered 0-99.
  When it is read into GRanges, 
  both the chromStart and chromEnd are on 1-based coordinate, \emph{i.e.},
  "chromStart"=1 and "chromEnd"=100.
}

\value{
  A \code{GRanges} is returned. 
  When no strand information is availabe in bed file,
  all the ranges are assumed to be on the positive strand.
}

\author{
  Ge Tan
}

\references{
  \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}
}

\seealso{
\code{\link{import.bed}}
}

\examples{
  bedHg19Fn <- file.path(system.file("extdata", package="CNEr"), 
                         "filter_regions.hg19.bed")
  bedHg19 <- readBed(bedHg19Fn)
}
