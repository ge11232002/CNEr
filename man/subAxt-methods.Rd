\name{subAxt-methods}
\docType{methods}
\alias{subAxt}
\alias{subAxt-methods}
\alias{subAxt,Axt-method}
\title{ \code{subAxt} method}
\description{
  Get subset of Axt alignments based on chromosome and ranges.
}

\usage{
  subAxt(x, chr, start, end,
      select=c("target", "query"), type=c("any", "within"),
      qSize=NULL)
}

\arguments{
  \item{x}{
    A object of \code{Axt}.
  }
  \item{chr}{
    A object of \code{character}. The chromosome name to extract.
  }
  \item{start, end}{
    A object of \code{integer}. 
    The range to extract is based on the positive strand.
  }
  \item{select}{
    When select is "target", 
    the subset criteria is for target alignments in axts.
    When select is "query", 
    the subset criteria is for query alignments in axts.
  }
  \item{type}{
    By default, any overlap is accepted. 
    If "type" is "within", 
    the selected axts must be wholly contained within the subject interval.
  }
  \item{qSize}{
    When select is "query", 
    qSize must be provided and is the length of chromosome chr.
  }
}

\details{
  Usually when we want to subset some axts from a \code{Axt} object,
  we care about all the axts within certain range.
  The axts can come from the axt file with chr as reference 
  (i.e., target sequence),
  or the axt file with chr as query sequence.
  When the chr is query sequence, it can be on the negative strand.
  Hence, the size of chromosome is necessary to 
  convert the search range to a range on negative strand coordinate.
}

\value{
  A subset of \code{Axt} object is returned.
}

\author{
  Ge Tan
}

\examples{
  axtFilesHg19DanRer7 = file.path(system.file("extdata", package="CNEr"), 
                                  "hg19.danRer7.net.axt")
  axtHg19DanRer7 = readAxt(axtFilesHg19DanRer7)
  subAxt(axtHg19DanRer7, chr="chr11", start=31500000, end=32500000,
         select="target", type="any")
}

\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
