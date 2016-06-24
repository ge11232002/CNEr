\name{plotCNEWidth}
\alias{plotCNEWidth}

\title{
  Plot the CNE width distribution.
}
\description{
  CNE widths can follow heavy tailed distribution.
  This function plot the reverse cumulative density distribution of CNE widths,
  and fit a discrete power-law distribution.
  Goodness of fit can also be evaluated.
}
\usage{
  plotCNEWidth(x, powerLawTest = FALSE, mc.cores=1L, ...)
}
  
\arguments{
  \item{x}{
    \code{GRangePairs} object: a pair of CNEs.
  }
  \item{powerLawTest}{
    \code{boolean}: when TRUE, the goodness of power law fit is evaluated.
  }
  \item{mc.cores}{
    \code{integer}: the number of threads to use in evalusation of 
    goodness of fit.
  }
  \item{\dots}{
    Additional points passed to \code{plot} function.
  }
}
\details{
  The power law distribution is usually for heavy tailed distribution.
  
  A reverse cumulative density distribution plot will be generated with optimal
  lower bound \emph{xmin}, scaling parameter\emph{alpha} for power law fit and 
  a pvalue of goodness of fit for power law distribution.
}
\value{
  An invisible list of fitted model is returned.
}
\references{
  Salerno, W., Havlak, P., and Miller, J. (2006). Scale-invariant structure of strongly conserved sequence in genomic intersections and alignments. Proc. Natl. Acad. Sci. U.S.A. 103, 13121-13125.
}
\author{
  Ge Tan
}
\note{
  The power law distribution implementation is based on package \pkg{poweRlaw}.
}

\examples{
  dbName <- file.path(system.file("extdata", package="CNEr"),
                      "danRer10CNE.sqlite")
  cneGRangePairs <- readCNERangesFromSQLite(dbName=dbName, 
                                            tableName="danRer10_hg38_45_50")
  plotCNEWidth(cneGRangePairs, powerLawTest=FALSE)
}