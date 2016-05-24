\name{readAxt}
\alias{readAxt}
\title{
  readAxt
}
\description{
  This function reads the \emph{axt} files into a \code{\link{Axt}} object.
}
\usage{
  readAxt(axtFiles)
}
\arguments{
  \item{axtFiles}{
    \code{character}(n): file names of the \emph{axt} files to read.
  }
}

\details{
  This function reads the \emph{axt} files of two assemblies.
  It can be a single big \emph{axt} file or several small \emph{axt} files.
  Different from the start coordinate in \emph{axt} file, 
  the start coordinate in \code{Axt} object is 1-based.
}

\value{
  A object \code{\link{Axt}} is returned.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{Axt}}
}

\examples{
  axtFile <- file.path(system.file("extdata", package="CNEr"), 
                       "hg38.danRer10.net.axt")
  axt <- readAxt(axtFile)
}
