\name{fetchChromSizes}
\alias{fetchChromSizes}
\title{
fetchChromSizes function.
}
\description{
This function tries to automate the fetch of chrom sizes for assembly 
from UCSC and other sources.
}
\usage{
fetchChromSizes(assembly)

}
\arguments{
  \item{assembly}{
  A \code{character} object: the canonical name of assembly, 
  i.e., hg19 for UCSC.
}
}
\details{
This function utilises mysql query for UCSC assemblies.
}
\value{
A object of \code{Seqinfo} is returned.
}
\author{
Ge Tan
}
\note{
Currently the assemblies from UCSC are supported.
}


\examples{
  fetchChromSizes("hg19")
  fetchChromSizes("mm10")
}

