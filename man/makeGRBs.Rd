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
           background=c("chromosome", "genome"), minCNEs=1L)
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
    can be "chromosome" or "genome". When using \code{slice} for the CNE density,
    the background is calculated on a per-chromosome or whole-genome basis.
  }
  \item{minCNEs}{
    \code{integer}(1): the minimal number of CNEs that a GRB needs to have.
  }
}
\details{
  First we calculate the CNE densities from the CNEs.
  Then we segment the regions according to the values of CNE densities.
  The regions with CNE densities above the expected CNE densities * ratio are
  considered as putative GRBs.
  Putative GRBs that do not encompass any gene are filtered out.
  Finally, the GRBs that have fewer than \code{minCNEs} number of CNEs will
  be filtered out.
}
\value{
  A \code{GRanges} object of GRB coordinates is returned.
  The numbers of CNEs and the coordinates of CNEs within each GRB are
  returned as a metadata column.
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
           background="chromosome", minCNEs=3L)
}
