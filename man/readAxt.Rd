\name{readAxt}
\alias{readAxt}
\title{
readAxt function.
}
\description{
This function reads the axt files into a \code{Axt} object.
}
\usage{
readAxt(axtFiles)
}
\arguments{
  \item{axtFiles}{
  Object of \code{character}. The length can be one or more.
}
}
\details{
  This function reads the axt files. 
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
  axtFilesHg19DanRer7 = file.path(system.file("extdata", package="CNEr"), 
                                  "hg19.danRer7.net.axt")
  axtHg19DanRer7 = readAxt(axtFilesHg19DanRer7)

}
