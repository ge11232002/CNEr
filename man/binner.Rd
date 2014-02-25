\name{binning-utils}
\alias{binFromCoordRange}
\alias{binRangesFromCoordRange}
\alias{binRestrictionString}

\title{UCSC bin indexing system utility functions}

\description{Utility functions for UCSC bin indexing system manipulation}

\usage{
  binFromCoordRange(starts, ends)
  binRangesFromCoordRange(start, end)
  binRestrictionString(start, end, field="bin")
}

\arguments{
  \item{starts, ends}{
    A vector of integers. A set of ranges.
  }
  \item{start, end}{
    A integer vector of length 1. A coordinate range.
  }
  \item{field}{
    Name of bin column. Default: "bin".
  }
}

\value{
  For \code{binFromCoordRange}, it returns the bin number 
  that should be assigned to a feature spanning the given range. 
  Usually it is used when creating a database for the features.

  For \code{binRangesFromCoordRange}, it returns the set of bin ranges 
  that overlap a given coordinate range. 
  It is usually used to find out the bins overlapped with a range.
  For SQL query, it is more convenient to use \code{binRestrictionString}
  than to use this function directly.

  For \code{binRestrictionString}, it returns a string to be used 
  in the WHERE section of a SQL SELECT statement 
  that is to select features overlapping a certain range. 
  * USE THIS WHEN QUERYING A DB *
}

\details{
  The UCSC bin indexing system was initially suggested by Richard Durbin and 
  Lincoln Stein to speed up the SELECT of a SQL query for the rows overlapping
  with certain genome coordinate.
  The system first used in UCSC genome browser is 
  described by Kent et. al. (2002).
}

\references{
Kent, W. J., Sugnet, C. W., Furey, T. S., Roskin, K. M., Pringle, 
T. H., Zahler, A. M., & Haussler, A. D. (2002). 
The Human Genome Browser at UCSC. Genome Research, 12(6), 
996-1006. doi:10.1101/gr.229102

\url{http://genomewiki.ucsc.edu/index.php/Bin_indexing_system}
}

\author{
  Ge Tan
}

\examples{
  binFromCoordRange(starts=c(10003, 1000000), ends=c(10004, 1100000))
  binRangesFromCoordRange(start=10000, end=2000000)
  binRestrictionString(start=10000, end=2000000, field="bin")
}
