\name{scoringMatrix}
\alias{scoringMatrix}

\title{
  scoringMatrix
}
\description{
  Generate the scoring matrix for \emph{lastz} aligner.
}
\usage{
  scoringMatrix(distance = c("far", "medium", "near"))
}

\arguments{
  \item{distance}{
    It can be "far", "medium" or "close". It decides the score matrix used in
    \emph{lastz} aligner.
    Generally, if two species are close to each other at human and chimp level,
    "close" should be used.
    If two species have a divergent time of 100 MYA, "far" should be used.
    In other cases, use "medium".
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
