\name{axtInfo}
\alias{axtInfo}
\title{
  axtInfo function
}
\description{
  Given the path of axt file, retrieve the alignments' withs information.
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
  It stores the withds of all the alignments.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{readAxt}}
}

\examples{
  axtFilesHg19DanRer7 <- file.path(system.file("extdata", package="CNEr"), 
                                   "hg19.danRer7.net.axt")
  axtInfoHg19DanRer7 <- axtInfo(axtFilesHg19DanRer7)
}
