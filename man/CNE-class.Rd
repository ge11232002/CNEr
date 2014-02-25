\name{CNE-class}
\Rdversion{1.1}
\docType{class}
\alias{CNE}
\alias{CNE-class}
\alias{assembly1}
\alias{assembly1,CNE-method}
\alias{assembly2}
\alias{assembly2,CNE-method}
\alias{CNE1}
\alias{CNE1,CNE-method}
\alias{CNE2}
\alias{CNE2,CNE-method}
\alias{CNEMerged}
\alias{CNEMerged,CNE-method}
\alias{CNERepeatsFiltered}
\alias{CNERepeatsFiltered,CNE-method}
\alias{thresholds}
\alias{thresholds,CNE-method}

\title{Class \code{"CNE"}}
\description{
  This class is used to store all intermediate and final results
  of CNE.
}

\usage{
### Constructors:
CNE(assembly1=character(), assembly2=character(),thresholds=character(),
    CNE1=list(), CNE2=list(), CNEMerged=list(), CNERepeatsFiltered=list(),
    alignMethod=character())

### Accessor-like methods:
\S4method{assembly1}{CNE}(x)
\S4method{assembly2}{CNE}(x)
\S4method{thresholds}{CNE}(x)
\S4method{CNE1}{CNE}(x)
\S4method{CNE2}{CNE}(x)
\S4method{CNEMerged}{CNE}(x)
\S4method{CNERepeatsFiltered}{CNE}(x)

## ... and more (see Methods)
}

\arguments{
    \item{assembly1}{Object of class \code{"character"}: 
      The name of assembly1.}
    \item{assembly2}{Object of class \code{"character"}:
      The name of assembly2.}
    \item{thresholds}{Object of class \code{"character"}:
      The thresholds of CNE scan: window size and identity score
      in the form of "49_50".}
    \item{CNE1}{Object of class \code{"list"}:
      The preliminary CNEs from axt file with assembly1 as reference.}
    \item{CNE2}{Object of class \code{"list"}:
      The preliminart CNEs from axt file with assembly2 as reference.}
    \item{CNEMerged}{Object of class \code{"list"}:
      The CNEs after merging CNE1 and CNE2.}
    \item{CNERepeatsFiltered}{Object of class \code{"list"}:
      The CNEs after being realigned back to reference genome,
      with blat in current implementation.}
    \item{alignMethod}{Object of class \code{"character"}:
      The method to realign CNEs back to reference genome.}
    \item{x}{Object of class \code{"CNE"}:
      A \code{"CNE"} object.}
}

\section{Methods}{
  \describe{
    \item{assembly1}{\code{signature(x = "CNE")}: 
      Get the assembly1 name.}
    \item{assembly2}{\code{signature(x = "CNE")}:
      Get the assembly2 name.}
    \item{CNE1}{\code{signature(x = "CNE")}:
      Get the CNE1 results.}
    \item{CNE2}{\code{signature(x = "CNE")}:
      Get the CNE2 results.} 
    \item{CNEMerged}{\code{signature(x = "CNE")}:
      Get the merged CNE results.}
    \item{CNERepeatsFiltered}{\code{signature(x = "CNE")}:
      Get the final CNE results.}
    \item{thresholds}{\code{signature(x = "CNE")}:
      Get the thresholds used for scanning CNEs.}
	 }
}

\author{
  Ge Tan
}

\examples{
showClass("CNE")
}

\keyword{classes}
