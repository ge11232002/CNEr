\name{dotplotAxt}
\alias{dotplotAxt}
\title{
  Syntenic dotplot
}
\description{
  Syntenic dotplot for \code{Axt} alignemnt object.
}
\usage{
  dotplotAxt(axt, targetSeqinfo, querySeqinfo,
             targetChrs = NULL, queryChrs = NULL,
             col = c("blue", "red"))
}
\arguments{
  \item{axt}{
    \code{Axt} object: the whole genome pairwise alignment of
    two compared species.
  }
  \item{targetSeqinfo,querySeqinfo}{
    \code{Seqinfo} object for both the target and query genomes.
  }
  \item{targetChrs,queryChrs}{
    \code{character}(n): the chromosomes to compare.
  }
  \item{col}{
    \code{character}(2): the colours for positive and negative strands.
  }
}
\details{
  This syntenic dotplot is a type of scatter plot.
  Each axis represents a concatenated selected chromosomes laid end-to-end, 
  and each dot in the scatter-plot represents a putative homologous match 
  between the two genomes.
  These dotplots are used for whole genome comparisons within the same genome or
  across two genomes from different taxa in order to identify synteny.
}
\value{
  A \code{ggplot} object.
}
\author{
  Ge Tan
}
\note{
  For highly fragmented assemblies, the synteny is invisible on the dotplot.
}
\examples{
  fn <- file.path(system.file("extdata", package="CNEr"),
                  "chr4.hg19.galGal3.net.axt.gz")
  axt <- readAxt(fn)
  targetSeqinfo <- fetchChromSizes("hg19")
  querySeqinfo <- fetchChromSizes("galGal3")
  targetChrs <- c("chr4")
  queryChrs <- c("chr4")
  dotplotAxt(axt, targetSeqinfo, querySeqinfo,
             targetChrs=targetChrs, queryChrs=queryChrs)
}