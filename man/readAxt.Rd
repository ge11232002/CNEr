\name{readAxt}
\alias{readAxt}
\title{
  readAxt
}
\description{
  Read the axt files into a \code{Axt} object.
}
\usage{
  readAxt(axtFiles)
}
\arguments{
  \item{axtFiles}{
    The \code{character}(n) file names of the 'axt' files to read.
  }
}

\details{
  This function reads the axt files of two assemblies.
  It can be single big axt file or 
  The coordinates in \code{Axt} object is 1-based.
}

\value{
  A object \code{Axt} is returned.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{Axt}}
}

\examples{
  axtFilesHg19DanRer7 <- file.path(system.file("extdata", package="CNEr"), 
                                   "hg19.danRer7.net.axt")
  axtHg19DanRer7 <- readAxt(axtFilesHg19DanRer7)
}
