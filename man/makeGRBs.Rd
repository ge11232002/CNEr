\name{makeGRBs}
\alias{makeGRBs}

\title{
  makeGRBs
}
\description{
  Make Genomic Regulatory Blocks (GRBs) boundaries prediction from a set of CNEs.
}
\usage{
  makeGRBs(x, winSize=NULL, genes=NULL, ratio=1,
           background=c("chromosome", "genome"))
}
  
\arguments{
  \item{x}{
    \code{GRangesList} object of a set of CNEs to use.
  }
  \item{winSize}{
    \code{integer}: the smoothing window size for CNE densities in kb.
    This value depends on the genome size of the reference genome.
    A larger genome requires bigger window size.
    For instance, 300kb is the appropriate window size for the human genome.
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
  \item{background}{
    \code{character}(1):
    can be "chromosome" or "genome". When slice the CNE density,
    the background is calculated on a per-chromosome or whole-genome basis.
  }
}
\details{
  First we calculated the CNE densities from the CNEs.
  Then we segment the regions according to the values of CNE densities.
  The regions with CNE densities above the expected CNE densities * ratio are
  considered as putative GRBs.
  As a last step, the putative GRBs that do not encompass any gene
  are filtered out.
}
\value{
  A \code{GRanges} object of GRB coordinates is returned, with the number of CNEs
  within each GRB.
}

\author{
  Ge Tan
}

\examples{
  library(TxDb.Drerio.UCSC.danRer10.refGene)
  refGenesDanRer10 <- genes(TxDb.Drerio.UCSC.danRer10.refGene)
  ancoraCNEsFns <- file.path(system.file("extdata", package="CNEr"), 
                             c("cne2wBf_cypCar1_danRer10_100_100",
                               "cne2wBf_cteIde1_danRer10_100_100",
                               "cne2wBf_AstMex102_danRer10_48_50"))
  cneList <- do.call(GRangesList, 
                     lapply(ancoraCNEsFns, readAncora, assembly="danRer10"))
  names(cneList) <- c("Common carp", "Grass carp", "Blind cave fish")
  seqlengths(cneList) <- seqlengths(TxDb.Drerio.UCSC.danRer10.refGene)[
                           names(seqlengths(cneList))]
  makeGRBs(cneList, winSize=200, genes=refGenesDanRer10, ratio=1.2,
           background="genome")
  makeGRBs(cneList, winSize=200, genes=refGenesDanRer10, ratio=1.2,
           background="chromosome")
}
