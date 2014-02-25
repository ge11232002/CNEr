\name{cneMerge}
\alias{cneMerge}

\title{
  CNE merge function
}
\description{
  Remove the CNEs which overlap on both genomes.
}
\usage{
cneMerge(cne1, cne2)
}

\arguments{
  \item{cne1, cne2}{
  A object of \code{data.frame}. The result from \code{ceScan}.
}
}


\value{
  A \code{data.frame} of CNEs is returned.
  In this table, the order of columns are consistent with cne1.
  For instance, if cne1 has the first three columns for zebrafish 
  and next three columns for human,
  in the merged table, the first three columns are 
  still the coordinates for zebrafish
  while the next three columns are coordinates for human.
}

\author{
  Ge Tan
}




\examples{
data(CNEHg19DanRer7)
data(CNEDanRer7Hg19)
cneMergedDanRer7Hg19 = mapply(cneMerge, CNEDanRer7Hg19, CNEHg19DanRer7, 
                              SIMPLIFY=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
