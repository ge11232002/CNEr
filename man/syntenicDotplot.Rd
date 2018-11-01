\name{syntenicDotplot-methods}
\docType{methods}
\alias{syntenicDotplot}
\alias{syntenicDotplot,Axt-method}
\alias{syntenicDotplot,GRangePairs-method}

\title{
  Syntenic dotplot
}
\description{
  Syntenic dotplot for \code{Axt} alignment object or \code{GRangePairs}.
}
\usage{
  syntenicDotplot(x, firstSeqlengths=NULL, secondSeqlengths=NULL,
                  firstChrs=NULL, secondChrs=NULL,
                  col=c("blue", "red"), type=c("line", "dot"))
}
\arguments{
  \item{x}{
    \code{Axt} object: the whole genome pairwise alignment of
    two species under comparison or
    \code{GRangePairs} object.
  }
  \item{firstSeqlengths,secondSeqlengths}{
    \code{integer}(n): seqlengths for both the first (target) and 
    second (query) genomes.
    When NULL, the seqlengths must exist in x.
  }
  \item{firstChrs,secondChrs}{
    \code{character}(n): the chromosomes to compare.
  }
  \item{col}{
    \code{character}(2): the colours for positive and negative strands.
  }
  \item{type}{
    \dQuote{line} or \dQuote{dot} plot type: When plotting massive number of
    ranges, \dQuote{dot} should be used. 
    Otherwise, \dQuote{line} should be used.
  }
}
\details{
  This syntenic dotplot is a type of scatter plot for \code{Axt} object,
  and line plot for \code{GRangePairs} object.
  In the case of possibly massive number of \code{Axt} alignments, 
  the line plots will make it invisible at a large genome scale.
  
  Each axis represents concatenated selected chromosomes laid end-to-end, 
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
  library(GenomeInfoDb)
  library(BSgenome.Ggallus.UCSC.galGal3)
  library(BSgenome.Hsapiens.UCSC.hg19)
  ## dotplot for Axt object
  fn <- file.path(system.file("extdata", package="CNEr"),
                  "chr4.hg19.galGal3.net.axt.gz")
  axt <- readAxt(fn)
  firstSeqlengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
  secondSeqlengths <- seqlengths(BSgenome.Ggallus.UCSC.galGal3)
  firstChrs <- c("chr4")
  secondChrs <- c("chr4")
  syntenicDotplot(axt, firstSeqlengths, secondSeqlengths,
                  firstChrs=firstChrs, secondChrs=secondChrs,
                  type="dot")
  
  ## dotplot for GRangePairs object
  data(grangesPairsForDotplot)
  syntenicDotplot(grangesPairsForDotplot, type="line")
}