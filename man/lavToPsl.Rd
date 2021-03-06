\name{lavToPsl}
\alias{lavToPsl}

\title{
  lavToPsl
}
\description{
  Wrapper function of \command{lavToPsl}: Convert blastz lav to psl format.
  This function doesn't work on Windows platform since Kent utilities only
  support Linux and Unix platforms.
}
\usage{
  lavToPsl(lavs, psls=sub("\\\.lav$", ".psl", lavs, ignore.case = TRUE),
           removeLav=TRUE, binary="lavToPsl")
}

\arguments{
  \item{lavs}{
    \code{character}(n): file names of input \emph{lav} files.
  }
  \item{psls}{
    code{character}(n): file names of output \emph{psl} files. 
    By default, in the same folder of input lav files with same names.
  }
  \item{removeLav}{
    \code{boolean}: When \code{TRUE}, the input \emph{lavs} files will be 
    removed after the conversion.
  }
  \item{binary}{
    \code{character}(1): the name/filename of the binary \command{lavToPsl}
    to call.
  }
}

\value{
  \code{character}(n): the file names of output \emph{psl} files.
}
\references{
  \url{http://hgdownload.cse.ucsc.edu/admin/exe/}
}
\author{
  Ge Tan
}

\seealso{
  \code{\link{lastz}}
}

\examples{
  \dontrun{
  ## This example doesn't run because it requires lav files from previous steps 
  ## and external Kent utilities.
    lavs <- tools::list_files_with_exts(
              dir="/Users/gtan/OneDrive/Project/CSC/CNEr/axt", exts="lav")
    lavToPsl(lavs, removeLav=FALSE, binary="lavToPsl")
  }
}

