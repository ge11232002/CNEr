\name{ceScan-methods}
\docType{methods}
\alias{ceScan}
\alias{ceScan,Axt-method}
\alias{ceScan,character-method}
\alias{ceScan,CNE-method}

\title{ceScan function}
\description{
  This is the main function for conserved noncoding elements (CNEs) 
  identification.
}

\usage{
  ceScan(x, tFilter=NULL, qFilter=NULL, qSizes, window=50L, identity=50L)
}

\arguments{
  \item{x}{
    \code{Axt} object or \code{character}(n) object of Axt filenames.
  }
  \item{tFilter}{
    \code{GRanges} object or NULL: regions to filter out for target assembly.
  }
  \item{qFilter}{
    \code{GRanges} object or NULL: regions to filter out for query assembly.
  }
  \item{qSizes}{
    \code{Seqinfo} object or NULL: it contains the \code{seqnames} and 
    \code{seqlengths} for query genome.
    This argument can be NULL when \code{qFilter} is NULL.
  }
  \item{window}{
    \code{integer}(n): the window size of scanning CNEs.
    By default, it is 50L.
  }
  \item{identity}{
    \code{integer}(n): the minimal identity score over the scanning window.
    By default, it is 50L. 
  }
}

\details{
  ceScan scan the axts alignmnets and identify the CNEs.
  ceScan can accept axts in \code{Axt} object and 
  filter in \code{GRanges} object,
  or directly the axt files and bed files.
  When the axt files and bed files are ready for computation, 
  it is recommended to use them directly rather than read them into R first.
  
  The details of algorithm will given in the vignette.
}

\value{
  A \code{list} of \code{data.frame} is returned. 
  Each element of the list is for one threshold.
}

\author{
  Ge Tan
}

\examples{
  axtFnHg38DanRer10 <- file.path(system.file("extdata", package="CNEr"), 
                                 "hg38.danRer10.net.axt")
  axtHg38DanRer10 <- readAxt(axtFnHg38DanRer10)
  axtFnDanRer10Hg38 <- file.path(system.file("extdata", package="CNEr"), 
                                 "danRer10.hg38.net.axt")
  axtDanRer10Hg38 <- readAxt(axtFnDanRer10Hg38)
  bedHg38Fn <- file.path(system.file("extdata", package="CNEr"), 
                        "filter_regions.hg38.bed")
  bedHg38 <- readBed(bedHg38Fn)
  bedDanRer10Fn <- file.path(system.file("extdata", package="CNEr"), 
                             "filter_regions.danRer10.bed")
  bedDanRer10 <- readBed(bedDanRer10Fn)
  qSizesHg38 <- fetchChromSizes("hg38")
  qSizesDanRer10 <- fetchChromSizes("danRer10")
  
  ## Axt object
  windows <- c(50L, 50L, 50L)
  identities <- c(45L, 48L, 49L)
  CNEHg38DanRer10 <- ceScan(x=axtHg38DanRer10, tFilter=bedHg38,
                            qFilter=bedDanRer10, qSizes=qSizesDanRer10,
                            window=windows, identity=identities)
  CNEDanRer10Hg38 <- ceScan(x=axtDanRer10Hg38, tFilter=bedDanRer10,
                            qFilter=bedHg38, qSizes=qSizesHg38,
                            window=windows, identity=identities)                  
                            
  ## Axt file names
  CNEHg38DanRer10 <- ceScan(x=axtFnHg38DanRer10, tFilter=bedHg38,
                            qFilter=bedDanRer10, qSizes=qSizesDanRer10,
                            window=windows, identity=identities)
}