\name{chainPreNet}
\alias{chainPreNet}

\title{
  chainPreNet
}
\description{
  Wrapper function of \command{chainPreNet}:
  Remove chains that don't have a chance of being netted.
  This function doesn't work on Windows platform since Kent utilities only
  support Linux and Unix platform.
}
\usage{
  chainPreNet(allChain, assemblyTarget, assemblyQuery,
              allPreChain=paste0(sub("\\\.2bit$", "", basename(assemblyTarget),
                                     ignore.case = TRUE), ".", 
                                 sub("\\\.2bit$", "", basename(assemblyQuery),
                                     ignore.case = TRUE), ".all.pre.chain"),
              removeAllChain=TRUE, binary="chainPreNet")
}

\arguments{
  \item{allChain}{
    \code{character}(1): file names of input \emph{allChain} file.
  }
  \item{assemblyTarget}{
    \code{character}(1): the file name of target assembly \emph{twoBit} file.
  }
  \item{assemblyQuery}{
    \code{character}(1): the file name of query assembly \emph{twoBit} file.
  }
  \item{allPreChain}{
    \code{character}(1): file names of merged \emph{allPreChain} file.
  }
  \item{removeAllChain}{
    \code{boolean}: When \code{TRUE}, the input \emph{allChain} file will be 
    removed after the conversion.
  }
  \item{binary}{
    \code{character}(1): the name/filename of the binary \command{chainPreNet}
    to call.
  }
}

\value{
  \code{character}(1): the file names of merged \emph{allPreChain} file.
}
\references{
  \url{http://hgdownload.cse.ucsc.edu/admin/exe/}
}
\author{
  Ge Tan
}

\seealso{
  \code{\link{chainMergeSort}}
}
\examples{
  \dontrun{
  ## This example doesn't run because it requires two bit files and external
  ## Kent utilities.
    allChain <- file.path("/Users/gtan/OneDrive/Project/CSC/CNEr/axt",
                          "danRer10.hg38.all.chain")
    assemblyTarget <- "/Users/gtan/OneDrive/Project/CSC/CNEr/2bit/danRer10.2bit"
    assemblyQuery <- "/Users/gtan/OneDrive/Project/CSC/CNEr/2bit/hg38.2bit"
    chainPreNet(allChain, assemblyTarget, assemblyQuery,
                allPreChain=file.path(
                              "/Users/gtan/OneDrive/Project/CSC/CNEr/axt",
                                      paste0(sub("\\\.2bit$", "", 
                                                 basename(assemblyTarget),
                                                 ignore.case = TRUE), ".", 
                                             sub("\\\.2bit$", "",
                                                 basename(assemblyQuery),
                                                 ignore.case = TRUE),
                                                 ".all.pre.chain")),
                removeAllChain=FALSE, binary="chainPreNet")
  }
}
