\name{makeCNEDensity}
\alias{makeCNEDensity}

\title{
  Make \sQuote{Bed}, \sQuote{bedGraph} and \sQuote{BigWig} files
}
\description{
  Make \sQuote{Bed}, \sQuote{bedGraph}, \sQuote{BigWig} files 
  from \code{GRangePairs} for
  display in other Genome Browser.
}
\usage{
  makeCNEDensity(x, outputDir = ".", 
                 genomeFirst = "first", genomeSecond = "second",
                 threshold = "50_50", 
                 windowSizeFirst = 300L, windowSizeSecond = 300L)
}
\arguments{
  \item{x}{
    \code{GRangePairs} object of CNEs.
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
    in format of "50_50".
  }
  \item{windowSizeFirst,windowSizeSecond}{
    \code{integer}(1): the smoothing window size for generating the CNE density
    in kb.
  }
}
\details{
  The CNE density is defined as the percentage of regions covered by CNEs
  within the smoothing window.
}
\value{
  The filenames of output \sQuote{Bed}, \sQuote{bedGraph} and
  \sQuote{BigWig} files.
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
  dbName <- file.path(system.file("extdata", package="CNEr"),
                      "danRer10CNE.sqlite")
  qAssemblyFn <- file.path(system.file("extdata",
                                       package="BSgenome.Hsapiens.UCSC.hg38"),
                           "single_sequences.2bit")
  tAssemblyFn <- file.path(system.file("extdata",
                             package="BSgenome.Drerio.UCSC.danRer10"),
                           "single_sequences.2bit")
  cneGRangePairs <- readCNERangesFromSQLite(dbName=dbName, 
                                            tableName="danRer10_hg38_45_50",
                                            tAssemblyFn=tAssemblyFn,
                                            qAssemblyFn=qAssemblyFn)
  makeCNEDensity(cneGRangePairs[1:1000])
}