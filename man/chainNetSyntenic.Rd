\name{chainNetSyntenic}
\alias{chainNetSyntenic}

\title{
  chainNetSyntenic
}
\description{
  Wrapper function of \emph{chainNetSyntenic}:
  Make alignment nets out of chains and add synteny info to net.
  This function doesn't work on Windows platform since Kent utilities only
  support Linux and Unix platform.
}
\usage{
  chainNetSyntenic(allPreChain, assemblyTarget, assemblyQuery,
                   netSyntenicFile=paste0(sub("\\\.2bit$", "",
                                              basename(assemblyTarget),
                                              ignore.case = TRUE), ".",
                                          sub("\\\.2bit$", "",
                                              basename(assemblyQuery),
                                              ignore.case = TRUE),
                                              ".noClass.net"),
                   binaryChainNet="chainNet", binaryNetSyntenic="netSyntenic")
}

\arguments{
  \item{allPreChain}{
    \code{character}(1): file names of input \emph{allPreChain} file.
  }
  \item{assemblyTarget}{
    \code{character}(1): the file name of target assembly \emph{twoBit} file.
  }
  \item{assemblyQuery}{
    \code{character}(1): the file name of query assembly \emph{twoBit} file.
  }
  \item{netSyntenicFile}{
    \code{character}(1): file names of output \emph{netSyntenicFile} file.
  }
  \item{binaryChainNet}{
    \code{character}(1): the name/filename of the binary \command{chainNet}
    to call.
  }
  \item{binaryNetSyntenic}{
    \code{character}(1): the name/filename of the binary \command{netSyntenic}
    to call.
  }
}
\details{
  Add classification information using the database tables:
  actually this step is not necessary in this pipeline
  according to http://blog.gmane.org/gmane.science.biology.ucscgenome.general/month=20130301.
  The class information will only be used for Genome Browser.
  Since it needs some specific modification of the table names for 
  certain species, we skip this step now.
  If this step is done, then the generated \emph{class.net} is 
  the gzipped net file that you see in UCSC Downloads area.
}
\value{
  \code{character}(1): the file names of generated \emph{net} file.
}
\references{
  \url{http://hgdownload.cse.ucsc.edu/admin/exe/}
}
\author{
  Ge Tan
}

\seealso{
  \code{\link{chainPreNet}}
}
\examples{
  \dontrun{
  ## This example doesn't run because it requires two bit files and external
  ## Kent utilities.
    allPreChain <- file.path("/Users/gtan/OneDrive/Project/CSC/CNEr/axt",
                             "danRer10.hg38.all.pre.chain")
    assemblyTarget <- "/Users/gtan/OneDrive/Project/CSC/CNEr/2bit/danRer10.2bit"
    assemblyQuery <- "/Users/gtan/OneDrive/Project/CSC/CNEr/2bit/hg38.2bit"
    chainNetSyntenic(allPreChain, assemblyTarget, assemblyQuery,
                     netSyntenicFile=file.path(
                       "/Users/gtan/OneDrive/Project/CSC/CNEr/axt",
                                               paste0(sub("\\\.2bit$", "",
                                                      basename(assemblyTarget),
                                                      ignore.case = TRUE), ".",
                                                      sub("\\\.2bit$", "",
                                                      basename(assemblyQuery),
                                                      ignore.case = TRUE),
                                              ".noClass.net")),
                     binaryChainNet="chainNet", binaryNetSyntenic="netSyntenic")
  }
}
