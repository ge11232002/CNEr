\name{chainMergeSort}
\alias{chainMergeSort}

\title{
  chainMergeSort
}
\description{
  Wrapper function of \command{chainMergeSort}:
  Combines sorted files into a larger sorted file.
  This function doesn't work on Windows platform since Kent utilities only
  support Linux and Unix platforms.
}
\usage{
  chainMergeSort(chains, assemblyTarget, assemblyQuery,
                 allChain=paste0(sub("\\\.2bit$", "", basename(assemblyTarget),
                                    ignore.case=TRUE), ".", 
                                 sub("\\\.2bit$", "", basename(assemblyQuery), 
                                     ignore.case=TRUE), ".all.chain"),
                 removeChains=TRUE, binary="chainMergeSort")
}

\arguments{
  \item{chains}{
    \code{character}(n): file names of input \emph{chains} files.
  }
  \item{assemblyTarget}{
    \code{character}(1): the file name of target assembly \emph{twoBit} file.
  }
  \item{assemblyQuery}{
    \code{character}(1): the file name of query assembly \emph{twoBit} file.
  }
  \item{allChain}{
    \code{character}(1): file names of merged \emph{allChain} file.
  }
  \item{removeChains}{
    \code{boolean}: When \code{TRUE}, the input \emph{chains} files will be 
    removed after the conversion.
  }
  \item{binary}{
    \code{character}(1): the name/filename of the binary \command{chainMergeSort}
    to call.
  }
}
\details{
  This \emph{allChain} file is what we get from UCSC download, \emph{e.g.}, 
  \bold{hg19.danRer7.all.chain.gz}.
}
\value{
  \code{character}(1): the file names of merged \emph{allChain} file.
}
\references{
  \url{http://hgdownload.cse.ucsc.edu/admin/exe/}
}
\author{
  Ge Tan
}

\seealso{
  \code{\link{axtChain}}
}
\examples{
  \dontrun{
  ## This example doesn't run because it requires two bit files and external
  ## Kent utilities.
    chains <- tools::list_files_with_exts(
                dir="/Users/gtan/OneDrive/Project/CSC/CNEr/axt", exts="chain")
    assemblyTarget <- "/Users/gtan/OneDrive/Project/CSC/CNEr/2bit/danRer10.2bit"
    assemblyQuery <- "/Users/gtan/OneDrive/Project/CSC/CNEr/2bit/hg38.2bit"
    chainMergeSort(chains, assemblyTarget, assemblyQuery,
      allChain=file.path("/Users/gtan/OneDrive/Project/CSC/CNEr/axt",
                         paste0(sub("\\\.2bit$", "", basename(assemblyTarget),
                                    ignore.case=TRUE), ".", 
                                sub("\\\.2bit$", "", basename(assemblyQuery), 
                                    ignore.case=TRUE), ".all.chain")),
      removeChains=FALSE, binary="chainMergeSort")
  }
}
