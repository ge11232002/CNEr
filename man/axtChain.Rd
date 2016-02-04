\name{axtChain}
\alias{axtChain}

\title{
  axtChain
}
\description{
  Wrapper function of \command{axtChain}: chain together psl alignments.
  If two matching alignments next to each other are close enough,
  they are joined into one fragment.
  This function doesn't work on Windows platform since Kent utilities only
  support Linux and Unix platform.
}
\usage{
  axtChain(psls, chains=sub("\\\.psl$", ".chain", psls, ignore.case=TRUE),
           assemblyTarget, assemblyQuery,
           distance=c("far", "medium", "far"),
           removePsl=TRUE, binary="axtChain")
}

\arguments{
  \item{psls}{
    \code{character}(n): file names of input \emph{psl} files.
  }
  \item{chains}{
    \code{character}(n): file names of output \emph{chain} files. 
    By default, in the same folder of input lav files with same names.
  }
  \item{assemblyTarget}{
    \code{character}(1): the file name of target assembly \emph{twoBit} file.
  }
  \item{assemblyQuery}{
    \code{character}(1): the file name of query assembly \emph{twoBit} file.
  }
  \item{distance}{
    It can be "far", "medium" or "close". It decides the score matrix used in
    \emph{lastz} aligner. 
    See `?scoringMatrix` for more details.
  }
  \item{removePsl}{
    \code{boolean}: When \code{TRUE}, the input \emph{psl} files will be 
    removed from the conversion.
  }
  \item{binary}{
    \code{character}(1): the name/filename of the binary \command{axtChain}
    to call.
  }
}

\value{
  \code{character}(n): the file names of output \emph{chain} files.
}
\references{
  \url{http://hgdownload.cse.ucsc.edu/admin/exe/}
}
\author{
  Ge Tan
}

\seealso{
  \code{\link{lavToPsl}}
}
\examples{
  \dontrun{
  ## This example doesn't run because it requires two bit files and external
  ## Kent utilities.
    psls <- tools::list_files_with_exts(
              dir="/Users/gtan/OneDrive/Project/CSC/CNEr/axt", exts="psl")
    assemblyTarget <- "/Users/gtan/OneDrive/Project/CSC/CNEr/2bit/danRer10.2bit"
    assemblyQuery <- "/Users/gtan/OneDrive/Project/CSC/CNEr/2bit/hg38.2bit"
    axtChain(psls, assemblyTarget=assemblyTarget,
             assemblyQuery=assemblyQuery, distance="far",
             removePsl=FALSE, binary="axtChain")
  }
}
