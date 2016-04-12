\name{mismatchSummary}
\alias{mismatchSummary}
\alias{mismatchSummary,Axt-method}

\title{
  Utility functions related to Axt alignment
}
\description{
  A collection of different functions used to deal with \code{Axt} object.
}
\usage{
  mismatchSummary(x, ...) ## mismatch number and proportion
}

\arguments{
  \item{x}{
    An \code{Axt} object
  }
  \item{...}{
    Currently not used.
  }
}
\details{
  'mismatchSummary': a numeric vector giving the numner of mismatches and 
    the proportion of mismatches.
}


\author{
  Ge Tan
}


\examples{
  axtFilesHg19DanRer7 <- file.path(system.file("extdata", package="CNEr"),
                                   "hg19.danRer7.net.axt")
  axtHg19DanRer7 <- readAxt(axtFilesHg19DanRer7)
  mismatchSummary(axtHg19DanRer7)
}