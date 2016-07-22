\name{plotCNEDistribution}
\alias{plotCNEDistribution}
\title{
  Plot the CNE genomic location
}
\description{
  Plot the CNE genomic location distribution.
  It gives a overview whether the CNEs form clusters.
}
\usage{
  plotCNEDistribution(x, chrs = NULL, chrScale = 1e+06)
}

\arguments{
  \item{x}{
    \code{GRanges} object: the CNE locations.
  }
  \item{chrs}{
    \code{character}(n): the chromosomes to show. By default,
    the largest 6 chromosomes/scaffolds are selected.
  }
  \item{chrScale}{
    \code{numeric}(1): the chromosome/scaffold scale in the plot. By default,
    it's \code{1e6} which is in \sQuote{Mb}.
  }
}
\details{
  In the plot, x axis is the genomic location along each chromosome/scaffold.
  The y axis is the sequential CNE number.
  A typical CNE cluster can be spotted by the dramatic increase in y axis and 
  small increase in x axis.
}
\value{
  A \code{ggplot} object.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{plotCNEWidth}}
}
\examples{
  library(rtracklayer)
  dbName <- file.path(system.file("extdata", package="CNEr"),
                      "danRer10CNE.sqlite")
  cneGRangePairs <- readCNERangesFromSQLite(dbName=dbName, 
                                            tableName="danRer10_hg38_45_50")
  twoBitFn <- file.path(system.file("extdata",
                                    package="BSgenome.Drerio.UCSC.danRer10"),
                        "single_sequences.2bit")
  x <- first(cneGRangePairs)
  seqlengths(x) <- seqlengths(TwoBitFile(twoBitFn))[names(seqlengths(x))]
  plotCNEDistribution(x) 
}