\name{summary}
\alias{summary}
\alias{summary,Axt-method}

\title{
  Utility functions related to Axt alignment
}
\description{
  A collection of different functions used to deal with \code{Axt} object.
}
\usage{
  summary(object, ...) ## mismatch number and proportion
}

\arguments{
  \item{object}{
    An \code{Axt} object
  }
  \item{...}{
    Currently not used.
  }
}
\value{
  A \code{table} object with the counts of mismatches, insertions, deletions
  and the matches of each base.
}
\author{
  Ge Tan
}

\examples{
  axtFilesHg38DanRer10 <- file.path(system.file("extdata", package="CNEr"),
                                   "hg38.danRer10.net.axt")
  axtHg38DanRer10 <- readAxt(axtFilesHg38DanRer10)
  summary(axtHg38DanRer10)
}
