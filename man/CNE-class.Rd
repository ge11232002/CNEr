\name{CNE-class}
\Rdversion{1.1}
\docType{class}
\alias{CNE}
\alias{CNE-class}
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
CNE(assembly1Fn=character(1), assembly2Fn=character(1),
    axt12Fn=character(), axt21Fn=character(),
    window=50L, identity=50L,
    CNE12=GRangePairs(), CNE21=GRangePairs(),
    CNEMerged=GRangePairs(), CNEFinal=GRangePairs(),
    aligner="blat", cutoffs1=4L, cutoffs2=4L)

### Accessor-like methods:
\S4method{thresholds}{CNE}(x)
\S4method{CNE12}{CNE}(x)
\S4method{CNE21}{CNE}(x)
\S4method{CNEMerged}{CNE}(x)
\S4method{CNEFinal}{CNE}(x)

## ... and more (see Methods)
}

\arguments{
    \item{assembly1Fn,assembly2Fn}{Object of class \code{"character"}: 
      The twoBit filenames of assembly1, assembly2}
    \item{axt12Fn,axt21Fn}{Object of class \code{"character"}: 
      The Axt filenames of assembly1 to assembly2, assembly2 to assembly1}
    \item{window}{Object of class \code{"integer"}:
      The window size for scanning CNEs. By default, it is 50.}
    \item{identity}{Object of class \code{"integer"}:
      The identity over the window size for scanning CNEs.
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
      The method to realign CNEs back to the reference genome.}
    \item{cutoffs1, cutoffs2}{Object of class \code{"integer"}:
      The CNEs with more than the cutoff hits on the reference genome 
      are removed.}
    \item{x}{Object of class \code{"CNE"}:
      A \code{"CNE"} object.}
}

\section{Methods}{
  \describe{
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
                       second=GRanges(seqnames=c("chr1"),
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
                       second=GRanges(seqnames=c("chr13", "chr4", "chr4"),
                                      ranges=IRanges(start=c(71727138,150679343,
                                                             146653164),
                                                     end=c(71727224, 150679400,
                                                           146653221)),
                                      strand="+")
                       )
  cne <- CNE(assembly1Fn=file.path(system.file("extdata",
                                    package="BSgenome.Drerio.UCSC.danRer10"),
                        "single_sequences.2bit"),
             assembly2Fn=file.path(system.file("extdata",
                                    package="BSgenome.Hsapiens.UCSC.hg38"),
                        "single_sequences.2bit"),
             window=50L, identity=50L,
             CNE12=CNE12, CNE21=CNE21, CNEMerged=CNE12, CNEFinal=CNE12,
             aligner="blat", cutoffs1=4L, cutoffs2=4L)
             
  ## Accessor
  CNE12(cne)
  CNE21(cne)
  thresholds(cne)
  CNEMerged(cne)
  CNEFinal(cne)
}
\keyword{classes}