\name{netToAxt}
\alias{netToAxt}

\title{
  netToAxt
}
\description{
  Wrapper function of \command{netToAxt} and \command{axtSort}:
  convert net (and chain) to axt, and sort axt files.
  This function doesn't work on the Windows platform since Kent utilities only
  support Linux and Unix platforms.
}
\usage{
  netToAxt(in.net, in.chain, assemblyTarget, assemblyQuery,
           axtFile=paste0(sub("\\\.2bit$", "", basename(assemblyTarget),
                              ignore.case = TRUE), ".",
                          sub("\\\.2bit$", "", basename(assemblyQuery),
                              ignore.case = TRUE), ".net.axt"),
           removeFiles=FALSE,
           binaryNetToAxt="netToAxt", binaryAxtSort="axtSort")
}

\arguments{
  \item{in.net}{
    \code{character}(1): file names of input \emph{net} file.
  }
  \item{in.chain}{
    \code{character}(1): file names of input \emph{chain} file.
  }
  \item{assemblyTarget}{
    \code{character}(1): the file name of target assembly \emph{twoBit} file.
  }
  \item{assemblyQuery}{
    \code{character}(1): the file name of query assembly \emph{twoBit} file.
  }
  \item{axtFile}{
    \code{character}(1): file names of output \emph{axt} file.
  }
  \item{removeFiles}{
    \code{boolean}: When \code{TRUE}, the input \emph{net} and \emph{chain}
    files will be removed after the conversion.
  }
  \item{binaryNetToAxt}{
    \code{character}(1): the name/filename of the binary \command{netToAxt}
    to call.
  }
  \item{binaryAxtSort}{
    \code{character}(1): the name/filename of the binary \command{axtSort}
    to call.
  }
}

\value{
  \code{character}(1): the file name of output \emph{axt} file.
}
\references{
  \url{http://hgdownload.cse.ucsc.edu/admin/exe/}
}
\author{
  Ge Tan
}

\seealso{
  \code{\link{chainNetSyntenic}}
}
\examples{
  \dontrun{
    ## This example doesn't run because it requires two bit files and external
    ## Kent utilities.
    in.net <- file.path("/Users/gtan/OneDrive/Project/CSC/CNEr/axt",
                        "danRer10.hg38.noClass.net")
    in.chain <- file.path("/Users/gtan/OneDrive/Project/CSC/CNEr/axt",
                          "danRer10.hg38.all.pre.chain")
    assemblyTarget <- "/Users/gtan/OneDrive/Project/CSC/CNEr/2bit/danRer10.2bit"
    assemblyQuery <- "/Users/gtan/OneDrive/Project/CSC/CNEr/2bit/hg38.2bit"
    netToAxt(in.net, in.chain, assemblyTarget, assemblyQuery,
             axtFile=file.path("/Users/gtan/OneDrive/Project/CSC/CNEr/axt",
                               paste0(sub("\\\.2bit$", "",
                                          basename(assemblyTarget),
                                          ignore.case = TRUE), ".",
                                      sub("\\\.2bit$", "",
                                          basename(assemblyQuery),
                                          ignore.case = TRUE),
                                      ".net.axt")),
             removeFiles=FALSE,
             binaryNetToAxt="netToAxt", binaryAxtSort="axtSort")
  }
}
