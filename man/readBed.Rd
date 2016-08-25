\name{readBed}
\alias{readBed}
\title{
  Read bed file
}
\description{
  Read the coordinates information from a bed file.
}
\usage{
  readBed(bedFile, assemblyFn=NULL)
}
\arguments{
  \item{bedFile}{
    \code{character}(1): filename of the \sQuote{bed} file to read.
  }
  \item{assemblyFn}{
    \code{character}(1): filename of the twoBit or fasta file of the genome.
  }
}
\details{
  This function is designed to read the bed file as
   \sQuote{chrom}, \sQuote{chromStart}, \sQuote{chromEnd}.
  The strand information is also stored where available.
  
  In the bed file, the \sQuote{chromStart} is on the 0-based coordinate system 
  while \sQuote{chromEnd} is on the 1-based coordinate system.
  For example, the first 100 bases of a chromosome are 
  defined as \sQuote{chromStart=0}, \sQuote{chromEnd=100}, 
  and span the bases numbered 0-99.
  When it is read into \code{GRanges},
  both the \sQuote{chromStart} and \sQuote{chromEnd} are on 1-based coordinate,
  \emph{i.e.},
  \sQuote{chromStart=1} and \sQuote{chromEnd=100}.
  
  When \sQuote{assemblyFn} is not \code{NULL}, the corresponding \code{Seqinfo}
  will be added into the returned \code{GRanges}.
}

\value{
  A \code{GRanges} object is returned. 
  When no strand information is available in the bed file,
  all the ranges are assumed to be on the positive strand.
}

\author{
  Ge Tan
}

\references{
  \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}
}

\seealso{
  \code{\link{import.bed}}
}

\examples{
  bedFn <- file.path(system.file("extdata", package="CNEr"), 
                     "filter_regions.hg38.bed")
  assemblyFn <- file.path(system.file("extdata",
                          package="BSgenome.Hsapiens.UCSC.hg38"),
                          "single_sequences.2bit")
  bed <- readBed(bedFn, assemblyFn=assemblyFn)
}
