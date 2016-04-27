\name{makeGRBs}
\alias{makeGRBs}

\title{
  makeGRBs
}
\description{
  Make Genomic Regulatory Blocks (GRBs) boundaries prediction from a set of CNEs.
}
\usage{
  makeGRBs(x, winSize=NULL, genes=NULL, ratio=0.5)
}
  
\arguments{
  \item{x}{
    \code{GRangesList} object of a set of CNEs to use.
  }
  \item{winSize}{
    \code{integer}: the smoothing window size for CNE densities in kb.
    This value depends on the genome size of the reference genome.
    A larger genome requires bigger windows size.
    For instance, 300kb is the appropriate windows size for human genome.
    By default, it is determined internally based on the genome size.
  }
  \item{genes}{
    \code{NULL} or \code{GRanges} object: the protein-coding genes ranges.
  }
  \item{ratio}{
    \code{numeric}(1): 
    the threshold to control the stringency of the GRBs.
    Higher value, shorter and fewer GRBs, and vice versa.
  }
}
\details{
  First we calculated the CNE densities from the CNEs.
  Then we segment the regions according to the value of CNE densities.
  The regions with CNE densities above the expected CNE densities * ratio are
  consider as putative GRBs.
  As last step, the putative GRBs that do not encompass any gene are filtered out.
}
\value{
  A \code{GRanges} object with GRB coordinates is returned.
}

\author{
  Ge Tan
}

\examples{
  \dontrun{
    ## Add example CNEs to make an example.
  }
}
