\name{blatCNE}
\alias{blatCNE}

\title{
Wrapper function of blat for CNEs
}
\description{
This wrapper function blat the CNEs against the reference genome.
}
\usage{
blatCNE(CNE, winSize, cutoffs1, cutoffs2, assembly1Twobit, assembly2Twobit,
        blatOptions=NULL, cutIdentity=90, tmpDir=tempdir(), blatBinary="blat")
}

\arguments{
  \item{CNE}{
    A object of \code{data.frame}. 
      Usually it is generated from \code{cneMerge} function.
  }
  \item{winSize}{
    A object of \code{integer}. 
    The window size used for identifying the CNEs, such as 50 or 30.
  }
  \item{cutoffs1, cutoffs2}{
    A object of \code{integer}.
    The CNEs with more than the cutoff hits on the reference genome are removed.
  }
  \item{assembly1Twobit, assembly2Twobit}{
    A object of \code{character}.
    The path of reference genome in two bit file format.
  }
  \item{blatOptions}{
    A object of \code{character}.
    When it is NULL, a bunch of preset parameters for blat 
    will be given based on the winSize parameter.
  }
  \item{cutIdentity}{
    A object of \code{integer}.
    Sets minimum sequence identity (in percent) in blat. Default is 90.
  }
  \item{tmpDir}{
    A object of \code{character}.
    By default, the R's temp dir is used.
    You can specify other path if your R's temp dir is small.
  }
  \item{blatBinary}{
    A object of \code{character}.
    The path of blat binary.
  }
}

\details{
  When winSize > 45, the blat options is 
  "-tileSize=11 -minScore=30 -repMatch=1024".
  
  When 35 < winSize <= 45, the blat options is 
  "-tileSize=10 -minScore=28 -repMatch=4096".
  
  When the winSize <= 35, the blat options is 
  "-tileSize=9 -minScore=24 -repMatch=16384".
}

\value{
  A \code{data.frame} containing the CNEs is returned.
}

\author{
Ge Tan
}

\examples{
  \dontrun{
    assemblyHg19Twobit = "/Users/gtan/CSC/CNEr/2bit/hg19.2bit"
    assemblyDanRer7Twobit = "/Users/gtan/CSC/CNEr/2bit/danRer7.2bit"
    cneBlatedDanRer7Hg19 = list()
    for(i in 1:length(cneMergedDanRer7Hg19)){
      cneBlatedDanRer7Hg19[[names(cneMergedDanRer7Hg19)[i]]] = 
      blatCNE(cneMergedDanRer7Hg19[[i]], 
        as.integer(sub("\\d+_", "", names(cneMergedDanRer7Hg19)[i])), 
        cutoffs1=4L, cutoffs2=8L, 
        assembly1Twobit=assemblyDanRer7Twobit, 
        assembly2Twobit=assemblyHg19Twobit,
        blatBinary="blat")
      }
    }
}


