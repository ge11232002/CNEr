\name{matchDistribution}
\alias{matchDistribution}
\alias{matchDistribution,Axt-method}
\title{
  Plot the distribution of matched alignments.
}
\description{
  Given a \code{Axt} alignment, plot a heatmap with percentage of
  each matched alignments.
}
\usage{
  matchDistribution(x, size=10000, title=NULL)
}

\arguments{
  \item{x}{
    \code{Axt} object.
  }
  \item{size}{
    \code{integer}(1): the number of alignments to use. By default, it is 10000.
  }
  \item{title}{
    \code{character}(1): the customised title for the plot.
  }
}
\details{
  By default, if there are more than 10,000 alignments, 
  for the consideration of speed,
  10,000 alignments will be sampled and calculated for the distribution.
  
  Only the four bases (A, C, G, T), gap (-) and any (N) are displayed. 
  Other ambiguous bases are not considere.
}
\value{
  A ggplot2 object will be returned.
}
\author{
  Ge Tan
}
\examples{
  axtFile <- file.path(system.file("extdata", package="CNEr"), 
                       "hg38.danRer10.net.axt")
  axt <- readAxt(axtFile)
  matchDistribution(axt)
}