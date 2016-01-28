\name{lavToPsl}
\alias{lavToPsl}

\title{
  lavToPsl
}
\description{
  Wrapper function of lavToPsl.
}
\usage{
  lavToPsl(lavs, psls = sub("\\.lav$", ".psl", lavs, ignore.case = TRUE),
           removeLav = TRUE)
}

\arguments{
  \item{lavs}{
    \code{character}(n): file names of input lav files
  }
  \item{psls}{
    code{character}(n): file names of output psl files. 
    By default, in the same folder of input lav files with same names.
  }
  \item{removeLav}{
    \code{boolean}: When \code{TRUE}, the input \code{lavs} files will be 
    removed from the conversion.
  }
}

\value{
  \code{character}(n): the file names of output psl files.
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
    lavs <- list_files_with_exts(dir="/Users/gtan/Downloads/2bit/axt",
                                 exts="lav")
    lavToPsl(lavs, removeLav=FALSE)
  }
}

