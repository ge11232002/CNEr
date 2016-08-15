\name{readAxt}
\alias{readAxt}
\title{
  Read \sQuote{Axt} file
}
\description{
  This function reads the \sQuote{Axt} files into a \code{\link{Axt}} object.
}
\usage{
  readAxt(axtFiles, tAssemblyFn=NULL, qAssemblyFn=NULL)
}
\arguments{
  \item{axtFiles}{
    \code{character}(n): filenames of the \sQuote{Axt} files to read.
  }
  \item{tAssemblyFn,qAssemblyFn}{
    \code{character}(1): filename of the \sQuote{twoBit} or \sQuote{fasta} 
    file for the target and query genome.
  }
}

\details{
  This function reads the \sQuote{Axt} files of two assemblies.
  It can be a single big \sQuote{Axt} file or several small \sQuote{Axt} files.
  Different from the start coordinate in \sQuote{Axt} file, 
  the start coordinate in \code{Axt} object is 1-based.
  
  When \sQuote{tAssemblyFn} and \sQuote{qAssemblyFn} are 
  not \code{NULL}, the corresponding \code{Seqinfo}
  will be added into the returned \code{Axt} object.
}

\value{
  A object \code{\link{Axt}} is returned.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{Axt}}
}

\examples{
  axtFile <- file.path(system.file("extdata", package="CNEr"), 
                       "hg38.danRer10.net.axt")
  tAssemblyFn <- file.path(system.file("extdata",
                           package="BSgenome.Hsapiens.UCSC.hg38"),
                           "single_sequences.2bit")
  qAssemblyFn <- file.path(system.file("extdata",
                           package="BSgenome.Drerio.UCSC.danRer10"),
                           "single_sequences.2bit")
  axt <- readAxt(axtFile, tAssemblyFn, qAssemblyFn)
}
