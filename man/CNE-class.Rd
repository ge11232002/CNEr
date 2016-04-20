\name{CNE-class}
\Rdversion{1.1}
\docType{class}
\alias{CNE}
\alias{CNE-class}
\alias{assembly1}
\alias{assembly1,CNE-method}
\alias{assembly2}
\alias{assembly2,CNE-method}
\alias{CNE12}
\alias{CNE12,CNE-method}
\alias{CNE21}
\alias{CNE21,CNE-method}
\alias{CNEMerged}
\alias{CNEMerged,CNE-method}
\alias{CNEFinal}
\alias{CNEFinal,CNE-method}
\alias{thresholds}
\alias{thresholds,CNE-method}

\title{Class \code{"CNE"}}
\description{
  \code{CNE} class contains all the meta-data of CNEs, including the pair of 
  assemblies, the thresholds, the intermediate and final CNE sets.
}

\usage{
### Constructors:
CNE(assembly1=character(1), assembly2=character(1),
    window=50L, identity=50L,
    CNE12=GRangePairs(), CNE21=GRangePairs(),
    CNEMerged=GRangePairs(), CNEFinal=GRangePairs(),
    aligner="blat", smoothingWindow1=300L, smoothingWindow2=300L)

### Accessor-like methods:
\S4method{assembly1}{CNE}(x)
\S4method{assembly2}{CNE}(x)
\S4method{thresholds}{CNE}(x)
\S4method{CNE12}{CNE}(x)
\S4method{CNE21}{CNE}(x)
\S4method{CNEMerged}{CNE}(x)
\S4method{CNEFinal}{CNE}(x)

## ... and more (see Methods)
}

\arguments{
    \item{assembly1}{Object of class \code{"character"}: 
      The name of assembly1.}
    \item{assembly2}{Object of class \code{"character"}:
      The name of assembly2.}
    \item{window}{Object of class \code{"integer"}:
      The window size of scanning CNEs. By default, it is 50.}
    \item{identity}{Object of class \code{"integer"}:
      The identity over the window size of scanning CNEs.
      By default, it is 50.}
    \item{CNE12}{Object of class \code{"GRangePairs"}:
      The preliminary CNEs from axt file with assembly1 as reference.}
    \item{CNE21}{Object of class \code{"GRangePairs"}:
      The preliminart CNEs from axt file with assembly2 as reference.}
    \item{CNEMerged}{Object of class \code{"GRangePairs"}:
      The CNEs after merging CNE1 and CNE2.}
    \item{CNEFinal}{Object of class \code{"GRangePairs"}:
      The CNEs after being realigned back to reference genome,
      with blat in current implementation.}
    \item{aligner}{Object of class \code{"character"}:
      The method to realign CNEs back to reference genome.}
    \item{x}{Object of class \code{"CNE"}:
      A \code{"CNE"} object.}
    \item{smoothingWindow1}{Object of class \code{"integer"}:
      The smoothing window for calculating the CNE density on assembly1.
      It is in kb and ranging from 10 to 1000.}
    \item{smoothingWindow2}{Object of class \code{"integer"}:
      The smoothing window for calculating the CNE density on assembly2.
      It is in kb and ranging from 10 to 1000.}
}

\section{Methods}{
  \describe{
    \item{assembly1}{\code{signature(x = "CNE")}: 
      Get the assembly1 name.}
    \item{assembly2}{\code{signature(x = "CNE")}:
      Get the assembly2 name.}
    \item{CNE12}{\code{signature(x = "CNE")}:
      Get the CNE1 results.}
    \item{CNE21}{\code{signature(x = "CNE")}:
      Get the CNE2 results.} 
    \item{CNEMerged}{\code{signature(x = "CNE")}:
      Get the merged CNE results.}
    \item{CNEFinal}{\code{signature(x = "CNE")}:
      Get the final CNE results.}
    \item{thresholds}{\code{signature(x = "CNE")}:
      Get the thresholds used for scanning CNEs.}
	 }
}

\author{
  Ge Tan
}

\examples{
  library(GenomicRanges)
  ## Constructor
  CNE12 <- GRangePairs(first=GRanges(seqnames=c("chr13", "chr4", "chr4"),
                                     ranges=IRanges(start=c(71727138,150679343,
                                                            146653164),
                                                    end=c(71727224, 150679400,
                                                          146653221)),
                                     strand="+"),
                       last=GRanges(seqnames=c("chr1"),
                                    ranges=IRanges(start=c(29854162, 23432387,
                                                           35711077),
                                                   end=c(29854248, 23432444,
                                                         35711134)),
                                    strand="+")
                       )
  CNE21 <- GRangePairs(first=GRanges(seqnames=c("chr1"),
                                    ranges=IRanges(start=c(29854162, 23432387,
                                                           35711077),
                                                   end=c(29854248, 23432444,
                                                         35711134)),
                                    strand="+"),
                       last=GRanges(seqnames=c("chr13", "chr4", "chr4"),
                                     ranges=IRanges(start=c(71727138,150679343,
                                                            146653164),
                                                    end=c(71727224, 150679400,
                                                          146653221)),
                                     strand="+")
                       )
  cne <- CNE(assembly1="hg38", assembly2="danRer10", window=50L, identity=50L,
             CNE12=CNE12, CNE21=CNE21, CNEMerged=CNE12, CNEFinal=CNE12,
             aligner="blat", smoothingWindow1=300L, smoothingWindow2=300L)
             
  ## Accessor
  assembly1(cne)
  assembly2(cne)
  CNE12(cne)
  CNE21(cne)
  thresholds(cne)
  CNEMerged(cne)
  CNEFinal(cne)
}

\keyword{classes}
