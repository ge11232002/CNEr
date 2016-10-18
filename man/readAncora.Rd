\name{readAncora}
\alias{readAncora}

\title{
  Read the cne file from Ancora format.
}
\description{
  Read the Ancora CNE file into a \code{GRanges} or \code{GRangePairs} object.
}
\usage{
  readAncora(fn, assembly=NULL, tAssemblyFn=NULL, qAssemblyFn=NULL)
}

\arguments{
  \item{fn}{
    \code{character}(1): the path of the Ancora CNE file in the format of 
    "cne2wBf_hg38_mm10_50_50".
  }
  \item{assembly}{
    \code{character}(1): the assembly to fetch. When it is NULL,
    CNEs on both assemblies are returned.
  }
  \item{tAssemblyFn,qAssemblyFn}{
    \code{character}(1): filename of the \sQuote{twoBit} or \sQuote{fasta} 
    file for the target and query genomes.
  }
}
\details{
  The Ancora CNE filename has its own naming style. For example, 
  "cne2wBf_hg38_mm10_50_50" denotes human coordinates for the first three
  columns of the file and mouse coordinates from the forth to the sixth column.
  
  The start coordinate system is 0-based.
}
\value{
  A \code{GRanges} object of the CNE ranges when assembly is specified,
  or a \code{GRangePairs} object when assembly is NULL.
}

\author{
  Ge Tan
}
\note{
  This function is mainly for internal use in Lenhard group.
}

\examples{
  fn <- file.path(system.file("extdata", package="CNEr"),
                  "cne2wBf_danRer10_hg38_45_50")
  zebrafishCNEs <- readAncora(fn, "danRer10")
  humanCNEs <- readAncora(fn, "hg38")
  zebrafishHumanCNEs <- readAncora(fn)
}