\name{saveCNEToSQLite-methods}
\docType{methods}
\alias{saveCNEToSQLite}
\alias{saveCNEToSQLite-methods}
\alias{saveCNEToSQLite,CNE,ANY,missing-method}
\alias{saveCNEToSQLite,data.frame,ANY,character-method}

\title{
saveCNEToSQLite function
}

\description{
This function save the CNE results into a local SQLite database.
}

\usage{
saveCNEToSQLite(CNE, dbName, tableName, overwrite=FALSE)
}

\arguments{
  \item{CNE}{
    An object of \code{data.frame}, the CNE data table
    or an object of \code{CNE}.
  }
  \item{dbName}{
    An object of \code{character}, the path of the local SQLite database.
  }
  \item{tableName}{
    An object of \code{character}, the name of table for this CNE data table,
    or missing when CNE is an object of \code{CNE}.
  }
  \item{overwrite}{
    An object of \code{boolean}, 
    whether or not to overwrite the table with same table name.
  }
}

\details{
  The input CNE table should have the colnames 
  "chr1", "start1", "end1", "chr2", "start2", "end2", 
  "strand", "similarity", "cigar".
  After the bin indexing, two additional columns "bin1" and "bin2" will be
  added before the column "chr1" and "chr2", respectively.
  
  If the input CNE is a \code{CNE} object,
  the tableName will be a combination of assembly names and thresholds.
  For instance, "danRer7_hg19_49_50" for "hg19" and "danRer7" with
  threshold "49_50".
}


\author{
Ge Tan
}

\examples{
  dbName = tempfile()
  data(cneBlatedDanRer7Hg19)
  for(i in 1:length(cneBlatedDanRer7Hg19)){
    tableName = paste("danRer7_hg19", names(cneBlatedDanRer7Hg19)[i],
                      sep="_")
    saveCNEToSQLite(cneBlatedDanRer7Hg19[[i]], dbName, tableName, 
                    overwrite=TRUE)
  }
  data(finalCNE)
  saveCNEToSQLite(finalCNE, dbName=dbName, overwrite=TRUE)
}

