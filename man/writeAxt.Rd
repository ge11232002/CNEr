\name{writeAxt}
\alias{writeAxt}

\title{
  \code{writeAxt} function
}
\description{
  Write an axt object into a file.
}
\usage{
  writeAxt(axt, con)
}

\arguments{
  \item{axt}{
    An \code{Axt} object to write.
  }
  \item{con}{
    A \link{connection} object or a character string.
  }
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
  axt <- readAxt(axtFile)
  writeAxt(axt, con=tempfile())
}



