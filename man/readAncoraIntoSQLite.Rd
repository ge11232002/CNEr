\name{readAncoraIntoSQLite}
\alias{readAncoraIntoSQLite}

\title{
  Read Ancora legacy CNE format
}
\description{
  Read Ancora legacy CNE format into a SQLite database.
}
\usage{
  readAncoraIntoSQLite(cneFns, dbName, overwrite=FALSE)
}

\arguments{
  \item{cneFns}{
    \code{character}(n): filenames of Ancora CNE files.
  }
  \item{dbName}{
    \code{character}(1): filename of SQLite database.
  }
  \item{overwrite}{
    \code{boolean}(1): whether or not to overwrite the existing table.
  }
}
\details{
  The Ancora legacy CNE file has the filename in the format of 
  "cne2wBf_AstMex102_danRer10_48_50".
  The first six columns are the coordinates of pairs of CNEs.
  The start coordinate system is 0-based and is converted into 1-based when
  it is imported into the SQLite database.
}
\value{
  A \code{character} vector of table names.
}
\author{
  Ge Tan
}
\note{
  This function is mainly for internal use in Lenhard group.
}

\seealso{
  \code{\link{readAncora}}
}
\examples{
  ancoraCNEsFns <- file.path(system.file("extdata", package="CNEr"),
                             c("cne2wBf_cypCar1_danRer10_100_100",
                               "cne2wBf_cteIde1_danRer10_100_100",
                               "cne2wBf_AstMex102_danRer10_48_50"))
  dbName <- tempfile()
  readAncoraIntoSQLite(ancoraCNEsFns, dbName, overwrite=FALSE)
}