\name{saveCNEToSQLite-methods}
\docType{methods}
\alias{saveCNEToSQLite}

\title{
  Save CNE to SQLite
}

\description{
  This function save the CNE results into a local SQLite database.
}

\usage{
  saveCNEToSQLite(x, dbName, tableName=NULL, overwrite=FALSE)
}

\arguments{
  \item{x}{
    An object of \code{CNE}, with CNEFinal computed or 
    a \code{GRangePairs} object.
  }
  \item{dbName}{
    \code{character}(1): the filename of the local SQLite database.
  }
  \item{tableName}{
    \code{character}(1): the name of table for this CNE data table.
    When it is NULL, the table name will be infered from the assembly filenames
    and scanning window/identity, in a format of "danRer10_hg38_49_50".
  }
  \item{overwrite}{
    \code{boolean}(1):
    whether or not to overwrite the existing table.
  }
}

\details{
  before loading into a SQLite database, 
  a bin indexing system is used to index the CNE range,
  which provides faster SQL query.
}


\author{
  Ge Tan
}

\examples{
  dbName <- tempfile()
  data(cneFinalListDanRer10Hg38)
  tableNames <- paste("danRer10", "hg38", names(cneFinalListDanRer10Hg38),
                      sep="_")
  for(i in 1:length(cneFinalListDanRer10Hg38)){
    saveCNEToSQLite(cneFinalListDanRer10Hg38[[i]], dbName, tableNames[i], 
                    overwrite=TRUE)
  }
}

