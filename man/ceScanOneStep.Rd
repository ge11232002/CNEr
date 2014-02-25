\name{ceScanOneStep}
\alias{ceScanOneStep}

\title{
ceScanOneStep function
}
\description{
This function run cne detection in one function.
}
\usage{
ceScanOneStep(axt1, filter1=NULL, sizes1, assembly1, twoBit1,
              axt2, filter2=NULL, sizes2, assembly2, twoBit2,
              thresholds=c("49_50"), blatBinary="blat",
              blatCutoff1, blatCutoff2)
}

\arguments{
  \item{axt1,axt2}{
  The \code{axt} object or axt filenames with each assembly as referecne.
 }
 \item{filter1,filter2}{
  The \code{GRanges} object or bed filenames.
 }
 \item{sizes1,sizes2}{
  A \code{Seqinfo} object which contains the \code{seqnames} and
    \code{seqlengths} for each assembly.
 }
 \item{assembly1,assembly2}{
  The assembly names.
 }
 \item{twoBit1,twoBit2}{
  The file names of two bit files of two assemblies.
 }
 \item{thresholds}{
    A \code{character} object specifiying the scanning windows and
    minimal score.
    It can be specified in th form of "45_50" with scanning windows 50
    and minial score 45.
    More than one thresholds can be provided.
 }
 \item{blatBinary}{
 A object of \code{character}. The path of blat binary.
 }
 \item{blatCutoff1, blatCutoff2}{
  A object of \code{integer}.
  The CNEs with more than the cutoff hits on the reference genome are removed.
 }
}


\value{
  An object \code{CNE} is returned.
}

\author{
  Ge Tan
}

