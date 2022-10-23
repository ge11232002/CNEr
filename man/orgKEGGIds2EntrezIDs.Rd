\name{orgKEGGIds2EntrezIDs}
\alias{orgKEGGIds2EntrezIDs}
\title{
  Fetch mapping from KEGG IDs to Entrez IDs
}
\description{
  Given the desired organism name, fetch the mapping between KEGG IDs and Entrez
  gene IDs.
}
\usage{
  orgKEGGIds2EntrezIDs(organism="Homo sapiens")
}

\arguments{
  \item{organism}{
    \code{character}(1): the name of organism to query. It has to be available
    at http://rest.kegg.jp/list/organism.
  }
}
\value{
  A \code{list} of Entrez gene IDs with KEGG IDs as names.
}
\author{
  Ge Tan
}

\examples{
  \dontrun{
    orgKEGGIds2EntrezIDs(organism="Homo sapiens")
  }
}
