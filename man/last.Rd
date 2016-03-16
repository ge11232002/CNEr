\name{lastal}
\alias{lastal}

\title{
  lastal wrapper
}
\description{
  Wrapper function of \command{lastal} to do the pairwise whole genome alignment.
  This function doesn't work on Windows platform.
}
\usage{
  lastal(db, queryFn,
         outputFn=sub("\\\.(fa|fasta)$", ".maf", 
                      paste(basename(db), basename(queryFn), sep = ","),
                      ignore.case = TRUE),
         distance=c("far", "medium", "close"), binary="lastal",
         mc.cores=getOption("mc.cores", 2L), echoCommand=FALSE)
}

\arguments{
  \item{db}{
    \code{character}(1): the file name of target assembly's lastal index. 
  }
  \item{queryFn}{
    \code{character}(1): the file name of query assembly \emph{fasta} file.
  }
  \item{outputFn}{
    \code{character}(1): the file name of the output \emph{maf} file.
  }
  \item{distance}{
    It can be "far", "medium" or "close". It decides the score matrix used in
    \emph{lastz} aligner. 
    See `?scoringMatrix` for more details.
  }
  \item{binary}{
    \code{character}(1): the name/filename of the binary \command{lastal} to call.
  }
  \item{mc.cores}{
    \code{integer}(1): the number of threads to use.
    By default, \code{getOption("mc.cores", 2L)}.
  }
  \item{echoCommand}{
    \code{boolean}(1): When \code{TRUE}, 
    only the command to run \command{lastal} is returned.
  }
}
\value{
  A \code{character}(1) vector of ouput \emph{maf} file names.
}
\references{
  \url{http://last.cbrc.jp/}
}
\author{
  Ge Tan
}
\note{
  \command{lastal} aligner must be installed on the machine to use this function.
}

\seealso{
  \code{\link{lastz}}
}
\examples{
  \dontrun{
    assemblyDir <- "/Users/gtan/OneDrive/Project/CSC/CNEr/2bit"
    ## Build the lastdb index
    system2(command="lastdb", args=c("-c", file.path(assemblyDir, "danRer10"),
                                     file.path(assemblyDir, "danRer10.fa")))

    ## Run lastal aligner
    lastal(db=file.path(assemblyDir, "danRer10"),
         queryFn=file.path(assemblyDir, "hg38.fa"),
         outputFn=file.path(axtDir, "danRer10.hg38.maf"),
         distance="far", binary="lastal", mc.cores=4L)

    ## maf to psl 
    psls <- file.path(axtDir, "danRer10.hg38.psl")
    system2(command="maf-convert", 
            args=c("psl", file.path(axtDir, "danRer10.hg38.maf"),
                   ">", psls))
  }
}

