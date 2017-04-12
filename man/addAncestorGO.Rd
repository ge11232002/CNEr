\name{addAncestorGO}
\alias{addAncestorGO}
\title{
  Add ancestor GO IDs
}
\description{
  Given a list of GO IDs, add the corresponding ancestor GO IDs.
}
\usage{
  addAncestorGO(go)
}
\arguments{
  \item{go}{
    A \code{list} of GO IDs. The elements of the list can be empty.
  }
}
\details{
  The ancestor GO IDs for each GO ID are added to the elements.
}
\value{
  A \code{list} of GO IDs with their ancestor GO IDs.
}
\author{
  Ge Tan
}
\note{
  This function is mainly designed for processing the gff annotation generated
  from interproscan, where for each gene, a set of GO IDs are assigned.
  However, for GO enrichment analysis, we need a list of mapping from genes to
  the GO IDs and their ancestor GO IDs as well.
}

\examples{
  library(GO.db)
  go <- list(c("GO:0005215", "GO:0006810", "GO:0016020"), "GO:0016579")
  addAncestorGO(go)
}
