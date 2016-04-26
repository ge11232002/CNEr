\name{readAncora}
\alias{readAncora}

\title{
  Read the cne file from Ancora format.
}
\description{
  Read the Ancora CNE file into a \code{GRanges}.
}
\usage{
  readAncora(fn, assembly)
}

\arguments{
  \item{fn}{
    \code{character}(1): the path of Ancora CNE file, in the format of 
    "cne2wBf_hg38_mm10_50_50".
  }
  \item{assembly}{
    \code{character}(1): the assembly to fetch.
  }
}
\details{
  The Ancora CNE filename has its own naming style. 
  "cne2wBf_hg38_mm10_50_50" will have human coordinates on the first three
  columns of the file and mouse coordinates on the forth to sixth column.
  
  The start coordinate is 0-based.
}
\value{
  \code{GRange} object of the CNE ranges will be returned.
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
}