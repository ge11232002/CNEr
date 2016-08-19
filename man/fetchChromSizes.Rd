\name{fetchChromSizes}
\alias{fetchChromSizes}
\title{
fetchChromSizes function.
}
\description{
This function tries to automate the fetch of chromosome sizes for assemblies 
from UCSC.
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
This function utilises mysql queries for UCSC assemblies.
}
\value{
A object of \code{Seqinfo} is returned.
}
\author{
  Ge Tan
}
\note{
  Currently, the assemblies from UCSC are supported.
}


\examples{
  \donttest{
    fetchChromSizes("hg19")
    fetchChromSizes("mm10")
  }
}

