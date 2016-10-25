\name{scoringMatrix}
\alias{scoringMatrix}

\title{
  scoringMatrix
}
\description{
  Generates the scoring matrix for \emph{lastz} aligner.
}
\usage{
  scoringMatrix(distance = c("far", "medium", "near"))
}

\arguments{
  \item{distance}{
    It can be "far", "medium" or "close". It defines the scoring matrix used in
    \emph{lastz} aligner.
    Generally, if two species are close to each other, for example human and chimp,
    "close" should be used.
    If two species have a divergence time of 100 MYA, "far" should be used.
    In other cases, "medium" should be used.
  }
}
\value{
  A \code{matrix} of the scoring matrix is returned.
}
\references{
  \url{http://genomewiki.ucsc.edu/index.php/Hg38_17-way_conservation_lastz_parameters}
}
\author{
  Ge Tan
}
\note{
  HOXD70 is medium. HoxD55 is far. human-chimp.v2 is close.
}

\seealso{
  \code{\link{lastz}}
}
\examples{
  scoringMatrix(distance="far")
}
