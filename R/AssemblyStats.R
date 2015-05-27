#########################################################################
# File Name: AssemblyStats.R
# Author: Ge Tan
# mail: gtan@me.com
# Created Time: Sun Nov 16 17:23:28 2014
#########################################################################

### -----------------------------------------------------------------
### N50 calculation for a assembly
### Exported!
NXX <- function(filepath, XX=50){
  if(grepl("\\.2bit$", filepath, ignore.case=TRUE)){
    lengths <- seqlengths(TwoBitFile(filepath))
  }else if(grepl("(\\.fa$|\\.fasta$)", filepath, ignore.case=TRUE)){
    lengths <- fasta.info(filepath)
  }else{
    stop("The suffix can only be .2bit, .fa, .fasta!")
  }
  lengths <- as.numeric(sort(lengths, decreasing=TRUE))
  index <- which(cumsum(lengths) / sum(lengths) >= XX/100)[1]
  return(lengths[index])
}

N50 <- function(filepath){
  return(NXX(filepath, XX=50))
}

N90 <- function(filepath){
  return(NXX(filepath, XX=90))
}


