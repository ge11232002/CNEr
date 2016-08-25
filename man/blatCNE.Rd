\name{blatCNE}
\alias{blatCNE}

\title{
  Wrapper function of blat for \code{CNE} object
}

\description{
  This wrapper function blat the CNEs against the reference genome. Note that blat must be installed on your system.
}

\usage{
  blatCNE(cne, blatOptions=NULL, cutIdentity=90)
}

\arguments{
  \item{cne}{
    \code{cne} object after \code{cneMerge} step.
  }
  \item{blatOptions}{
    \code{character}(1): the \command{blat} options.
    When it is NULL, the options will be chosen based on the window size for
    scanning CNEs.
  }
  \item{cutIdentity}{
    \code{integer}(1): the minimum sequence identity (in percent) for a match
    in \command{blat}. By default, it is 90.
  }
}

\details{
  When winSize > 45, the blat option is 
  "-tileSize=11 -minScore=30 -repMatch=1024".
  
  When 35 < winSize <= 45, the blat option is 
  "-tileSize=10 -minScore=28 -repMatch=4096".
  
  When the winSize <= 35, the blat option is 
  "-tileSize=9 -minScore=24 -repMatch=16384".
}

\value{
  A \code{CNE} object with a final set of CNEs.
}

\author{
  Ge Tan
}

\examples{
  \dontrun{
    data(CNEDanRer10Hg38)
    data(CNEHg38DanRer10)
    cne <- CNE(assembly1Fn=file.path(system.file("extdata",
                                    package="BSgenome.Drerio.UCSC.danRer10"),
                        "single_sequences.2bit"),
               assembly2Fn=file.path(system.file("extdata",
                                    package="BSgenome.Hsapiens.UCSC.hg38"),
                        "single_sequences.2bit"),
               window=50L, identity=45L, CNE12=CNEDanRer10Hg38[["45_50"]],
               CNE21=CNEHg38DanRer10[["45_50"]], aligner="blat")
    cne <- cneMerge(cne)
    cne <- blatCNE(cne)
  }
}