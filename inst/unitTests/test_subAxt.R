

test_subAxt <- function(){
  data(axtHg19DanRer7)
  ## Check the nr of axts
  checkIdentical(length(axtHg19DanRer7), 133L)

  ## Check the subAxt on target
  foo <- subAxt(axtHg19DanRer7, chr="chr11", start=31500113, end=31500120, 
                select="target")
  checkIdentical(as.character(targetSeqs(foo)), "AAATGCAG")
  checkIdentical(as.character(querySeqs(foo)), "GAGTGC-T")

  ## Check the subAxt on multiple ranges
  foo <- subAxt(axtHg19DanRer7, chr="chr11", start=c(31082021, 32461267),
                end=c(31082862,32461581), select="target")
  checkIdentical(length(foo), 3L)
  
}

