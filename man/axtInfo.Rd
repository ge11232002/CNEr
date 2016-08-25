\name{axtInfo}
\alias{axtInfo}
\title{
  axtInfo function
}
\description{
  Given the path of the axt file, this function retrieves information on the widths of the alignments. 
}
\usage{
  axtInfo(axtFiles)
}

\arguments{
  \item{axtFiles}{
    The filenames of axt files.
  }
}

\value{
  A vector of \code{integer} is returned. 
  It stores the widths of all the alignments.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{readAxt}}
}

\examples{
  axtFile <- file.path(system.file("extdata", package="CNEr"), 
                       "hg38.danRer10.net.axt")
  axtInfo <- axtInfo(axtFile)
}
