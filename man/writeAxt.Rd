\name{writeAxt}
\alias{writeAxt}

\title{
  \code{writeAxt} function
}
\description{
  Write an axt object into file.
}
\usage{
writeAxt(axt, con)
}

\arguments{
  \item{axt}{
    A \code{Axt} object to be written.
  }
  \item{con}{
    A connection object or a character string.
  }
}


\author{
  Ge Tan
}

\seealso{
  \code{\link{readAxt}}
}
\examples{
  axtFilesHg19DanRer7 = file.path(system.file("extdata", package="CNEr"),
                                  "hg19.danRer7.net.axt")
  axtHg19DanRer7 = readAxt(axtFilesHg19DanRer7)
  writeAxt(axtHg19DanRer7, con=tempfile())
}



