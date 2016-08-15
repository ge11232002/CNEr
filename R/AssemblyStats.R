### -----------------------------------------------------------------
### N50 calculation for a assembly
### Exported!
NXX <- function(filepath, XX=50){
  if(file_ext(filepath) == "2bit"){
    lengths <- seqlengths(TwoBitFile(filepath))
  }else if(file_ext(filepath) %in% c("fa", "fasta")){
    lengths <- fasta.seqlengths(filepath)
  }else{
    stop("The suffix can only be .2bit, .fa, .fasta!")
  }
  lengths <- as.numeric(sort(lengths, decreasing=TRUE))
  index <- which(cumsum(lengths) / sum(lengths) >= XX/100)[1]
  return(lengths[index])
}

N50 <- function(fn){
  return(NXX(fn, XX=50))
}

N90 <- function(fn){
  return(NXX(fn, XX=90))
}


