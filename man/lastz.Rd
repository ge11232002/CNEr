\name{lastz}
\alias{lastz}

\title{
  lastz wrapper
}
\description{
  Wrapper function of lastz to do the pairwise whole genome alignment.
}
\usage{
  lastz(assemblyTarget, assemblyQuery, outputDir = ".",
        chrsTarget = NULL, chrsQuery = NULL,
        distance = c("far", "medium", "close"), binary = "lastz", 
        mc.cores = getOption("mc.cores", 2L), echoCommand = FALSE)
}

\arguments{
  \item{assemblyTarget}{
    \code{character}(1): the file name of target assembly \emph{twoBit} file.
}
  \item{assemblyQuery}{
    \code{character}(1): the file name of query assembly \emph{twoBit} file.
}
  \item{outputDir}{
    \code{character}(1): the folder to put the generated \emph{lav} files.
}
  \item{chrsTarget}{
    NULL or \code{character}(n): when it's NULL, all the available chromosomes
    from the target assembly will be aligned.
}
  \item{chrsQuery}{
    NULL or \code{character}(n): when it's NULL, all the available chromosomes
    from the query assembly will be aligned.
}
  \item{distance}{
    It can be "far", "medium" or "close". It decides the score matrix used in
    \emph{lastz} aligner.
    Generally, if two species are close to each other at human and chimp level,
    "close" should be used.
    If two species have a divergent time of 100 MYA, "far" should be used.
    In other cases, use "medium".
}
  \item{binary}{
    \code{character}(1): the name/filename of the binary \emph{lastz} to call.
}
  \item{mc.cores}{
    \code{integer}(1): the number of threads to use. By default, \code{getOption("mc.cores", 2L)}.
}
  \item{echoCommand}{
    \code{boolean}(1): When \code{TRUE}, 
    only the command to run \emph{lastz} is returned.
  }
}

\value{
  A \code{character}(n) vector of ouput \emph{lav} file names.
}
\references{
  \url{http://www.bx.psu.edu/~rsharris/lastz/}
}
\author{
  Ge Tan
}
\note{
  \emph{lastz} aligner must be installed on the machine to use this function.
}

\seealso{
  \code{\link{lavToPsl}}
}

\examples{
  \dontrun{
    assemblyTarget <- "/Users/gtan/Downloads/2bit/fr3.2bit"
    assemblyQuery <- "/Users/gtan/Downloads/2bit/danRer10.2bit"
    lavs <- lastz(assemblyTarget, assemblyQuery, 
                  outputDir="/Users/gtan/Downloads/2bit/axt",
			            chrsTarget=c("chr1", "chr2", "chr3"),
			            chrsQuery=c("chr1", "chr2", "chr3"),
			            distance="far", mc.cores=4)
  }
}
