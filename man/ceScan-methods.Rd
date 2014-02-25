\name{ceScan-methods}
\docType{methods}
\alias{ceScan-methods}
\alias{ceScan}
\alias{ceScan,Axt,GRanges,GRanges,Seqinfo-method}
\alias{ceScan,Axt,GRanges,missing,missing-method}
\alias{ceScan,Axt,missing,GRanges,Seqinfo-method}
\alias{ceScan,Axt,missing,missing,missing-method}
\alias{ceScan,character,character,character,Seqinfo-method}
\alias{ceScan,character,character,missing,missing-method}
\alias{ceScan,character,missing,character,Seqinfo-method}
\alias{ceScan,character,missing,missing,missing-method}
\title{ceScan function}
\description{
  This is the main function for conserved noncoding elements (CNEs) 
  identification.
}

\usage{
  ceScan(axts, tFilter, qFilter, qSizes, thresholds="49_50")
}

\arguments{
  \item{axts}{
    A \code{Axt} object or \code{character} object with the paths of axt files.
  }
  \item{tFilter}{
    A \code{GRanges} object or \code{character} 
    object with the path of bed file for target genome filter. 
    This argument can also be missing when target filter is not available.
  }
  \item{qFilter}{
    A \code{GRanges} object or \code{character} 
    object with the path of bed file for
    query genome filter.
    This argument can also be missing when query filter is not available.
  }
  \item{qSizes}{
    A \code{Seqinfo} object which contains the \code{seqnames} and 
    \code{seqlengths} for query genome.
    This argument can be missing when \code{qFilter} is missing.
  }
  \item{thresholds}{
    A \code{character} object specifiying the scanning windows and 
    minimal score.
    It can be specified in th form of "45_50" with scanning windows 50 
    and minial score 45.
    More than one thresholds can be provided.
  }
}
\section{Methods}{
\describe{

\item{\code{signature(axts = "Axt", tFilter = "GRanges", qFilter = "GRanges", 
qSizes = "Seqinfo")}}{
}

\item{\code{signature(axts = "Axt", tFilter = "GRanges", qFilter = "missing", 
qSizes = "missing")}}{
}

\item{\code{signature(axts = "Axt", tFilter = "missing", qFilter = "GRanges", 
  qSizes = "Seqinfo")}}{
}

\item{\code{signature(axts = "Axt", tFilter = "missing", qFilter = "missing", 
  qSizes = "missing")}}{
}

\item{\code{signature(axts = "character", tFilter = "character", 
  qFilter = "character", qSizes = "Seqinfo")}}{
}

\item{\code{signature(axts = "character", tFilter = "character", 
  qFilter = "missing", qSizes = "missing")}}{
}

\item{\code{signature(axts = "character", tFilter = "missing", 
  qFilter = "character", qSizes = "Seqinfo")}}{
}

\item{\code{signature(axts = "character", tFilter = "missing", 
  qFilter = "missing", qSizes = "missing")}}{
}
}}

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
  axtFilesHg19DanRer7 = file.path(system.file("extdata", package="CNEr"), 
                                  "hg19.danRer7.net.axt")
  axtHg19DanRer7 = readAxt(axtFilesHg19DanRer7)
  axtFilesDanRer7Hg19 = file.path(system.file("extdata", package="CNEr"), 
                                  "danRer7.hg19.net.axt")
  axtDanRer7Hg19 = readAxt(axtFilesDanRer7Hg19)
  bedHg19Fn = file.path(system.file("extdata", package="CNEr"), 
                        "filter_regions.hg19.bed")
  bedHg19 = readBed(bedHg19Fn)
  bedDanRer7Fn = file.path(system.file("extdata", package="CNEr"), 
                           "filter_regions.danRer7.bed")
  bedDanRer7 = readBed(bedDanRer7Fn)
  qSizesHg19 = fetchChromSizes("hg19")
  qSizesDanRer7 = fetchChromSizes("danRer7")
  CNEHg19DanRer7 = ceScan(axts=axtHg19DanRer7, tFilter=bedHg19,
                          qFilter=bedDanRer7, qSizes=qSizesDanRer7,
                          thresholds=c("45_50", "48_50", "49_50"))
  CNEDanRer7Hg19 = ceScan(axts=axtDanRer7Hg19, tFilter=bedDanRer7,
                          qFilter=bedHg19, qSizes=qSizesHg19,
                          thresholds=c("45_50", "48_50", "49_50"))

}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
