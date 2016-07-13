\name{syntenicDotplot-methods}
\docType{methods}
\alias{syntenicDotplot}
\alias{syntenicDotplot,Axt-method}
\alias{syntenicDotplot,GRangePairs-method}

\title{
  Syntenic dotplot
}
\description{
  Syntenic dotplot for \code{Axt} alignemnt object or \code{GRangePairs}.
}
\usage{
  syntenicDotplot(x, firstSeqlengths=NULL, lastSeqlengths=NULL,
                  firstChrs=NULL, lastChrs=NULL,
                  col=c("blue", "red"))
}
\arguments{
  \item{x}{
    \code{Axt} object: the whole genome pairwise alignment of
    two compared species or
    \code{GRangePairs} object.
  }
  \item{firstSeqlengths,lastSeqlengths}{
    \code{integer}(n): seqlengths for both the first(target) and 
    last(query) genomes.
    When NULL, the seqlengths must exist in x.
  }
  \item{firstChrs,lastChrs}{
    \code{character}(n): the chromosomes to compare.
  }
  \item{col}{
    \code{character}(2): the colours for positive and negative strands.
  }
}
\details{
  This syntenic dotplot is a type of scatter plot for \code{Axt} object,
  and line plot for \code{GRangePairs} object.
  Because the possible massive number of \code{Axt} alignments, 
  the line plots will make it invisible on large genome scale.
  
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
  library(GenomeInfoDb)
  
  ## dotplot for Axt object
  fn <- file.path(system.file("extdata", package="CNEr"),
                  "chr4.hg19.galGal3.net.axt.gz")
  axt <- readAxt(fn)
  firstSeqlengths <- seqlengths(fetchChromSizes("hg19"))
  lastSeqlengths <- seqlengths(fetchChromSizes("galGal3"))
  firstChrs <- c("chr4")
  lastChrs <- c("chr4")
  syntenicDotplot(axt, firstSeqlengths, lastSeqlengths,
                  firstChrs=firstChrs, lastChrs=lastChrs)
  
  ## dotplot for GRangePairs object
  data(grangesPairsForDotplot)
  syntenicDotplot(grangesPairsForDotplot)
}