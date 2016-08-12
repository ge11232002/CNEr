\name{makeAncoraFiles}
\alias{makeAncoraFiles}

\title{
  makeAncoraFiles
}
\description{
  Make ancora format files from \code{GRangePairs} of CNE
}
\usage{
  makeAncoraFiles(cne, outputDir = ".", 
                  genomeFirst = "first", genomeSecond = "second",
                  threshold = "50_50")
}

\arguments{
  \item{cne}{
    \code{GRangePairs} object of CNE.
  }
  \item{outputDir}{
    \code{character}(1): the output directory of 
    \sQuote{Bed} and \sQuote{BigWig} files.
  }
  \item{genomeFirst,genomeSecond}{
    \code{character}(1): the genome name of the first and second species.
  }
  \item{threshold}{
    \code{character}(1): the threshold used to identify the CNEs
    in format of "50_50"..
  }
}
\value{
  The filenames of output.
}

\author{
  Ge Tan
}
\note{
  This function is mainly for internal use in Lenhard group.
}

\seealso{
  \code{\link{readAncora}}
}
\examples{
  data(cneFinalListDanRer10Hg38)
  cne <- CNEFinal(cneFinalListDanRer10Hg38[["45_50"]])
  makeAncoraFiles(cne, genomeFirst = "danRer10", genomeSecond = "hg38",
                  threshold = "45_50")
}