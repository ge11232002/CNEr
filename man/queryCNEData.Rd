\name{queryCNEData}
\alias{queryCNEData}

\title{
  Query the CNEData package to fetch the CNEs
}
\description{
  Query the CNEData package to fetch the CNEs based on target, query species,
  winSize and identity.
}
\usage{
  queryCNEData(dbName, target, query, winSize, identity,
               type=c("target", "all"))
}

\arguments{
  \item{dbName}{
    The path of SQLite database.
  }
  \item{target, query}{
    The CNEs between target and query species.
  }
  \item{winSize, identity}{
    The thresholds of CNEs to fetch on identity over winSize.
  }
  \item{type}{
    Which set of CNEs are returned. When it is "all", 
    the CNEs of target always on the left side of returned data.frame.
  }
}
\value{
  A data.frame of CNEs coordinates in chr, start, end.
}

\author{
  Ge Tan
}


\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

