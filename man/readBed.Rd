\name{readBed}
\alias{readBed}
\title{
readBed function
}
\description{
readBed reads bed file with an embeded C IO function.
}
\usage{
readBed(bedFile)
}
\arguments{
  \item{bedFile}{
  The path of bed file.
}
}
\details{
  This function is designed to read the bed file only 
  with three mandatory columns, 
  i.e., "chrom", "chromStart", "chromEnd".
  This function utilises the C interface for speedy import.
  For general bed file import, 
  please use the \code{import.bed} from package \code{rtracklayer}.
  
  In bed file, the "chromStart" is on 0-based coordinate 
  while "chromEnd" is on 1-based coordinate.
  For example, the first 100 bases of a chromosome are 
  defined as chromStart=0, chromEnd=100, 
  and span the bases numbered 0-99.
  When it is read into GRanges, 
  both the chromStart and chromEnd are on 1-based coordinate, i.e.,
  chromStart=1 and chromEnd=100.
}
\value{
  A \code{GRanges} is returned. 
  When no strand information is availabe in bed file,
  all the ranges are assumed to be on the positive strand.
}


\author{
Ge Tan
}


\seealso{
\code{\link{import.bed}}
}

\examples{
  bedHg19Fn = file.path(system.file("extdata", package="CNEr"), 
                        "filter_regions.hg19.bed")
  bedHg19 = readBed(bedHg19Fn)
}
