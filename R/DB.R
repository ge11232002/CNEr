### -----------------------------------------------------------------
### query the new CNEData data package
### Not exported!
#dbName = "/mnt/biggley/projects/gtan/projects/CNEr/CNErData/CNEData.sqlite"
#target = "equCab2"
#query = "canFam3"
#type = "all"
#winSize = 50
#identity = 49
#foo = queryCNEData(dbName, target, query, winSize, identity, type)

queryCNEData <- function(dbName, target, query, winSize, identity,
                         type=c("target", "all")
                         ){
  type <- match.arg(type)
  con <- dbConnect(SQLite(), dbname=dbName)
  on.exit(dbDisconnect(con))
  tableName <- ifelse(target < query, paste(target, query, sep="_"),
                      paste(query, target, sep="_"))
  tableName <- paste(tableName, identity, winSize, sep="_")
  if(type == "target"){
    if(target < query){
      sqlCmd <- "SELECT chr1, start1, end1"
    }else{
      sqlCmd <- "SELECT chr2, start2, end2"
    }
  #}else if(type == "query"){
  #  if(target < query){
  #    sqlCmd <- "SELECT chr2, start2, end2"
  #  }else{
  #    sqlCmd <- "SELECT chr1, start1, end1"
  #  }
  }else if(type == "all"){
    if(target < query){
      sqlCmd <- "SELECT chr1, start1, end1, chr2, start2, end2"
    }else{
      sqlCmd <- "SELECT chr2, start2, end2, chr1, start1, end1"
    }
  }
  sqlCmd <- paste(sqlCmd, "FROM", tableName)
  cne <- dbGetQuery(con, sqlCmd)
  return(cne)
}

