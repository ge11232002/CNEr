

test_subAxt <- function(){
  axtFilesHg19DanRer7 <- file.path(system.file("extdata", package="CNEr"),
                                   "hg19.danRer7.net.axt")
  axtHg19DanRer7 <- readAxt(axtFilesHg19DanRer7)
  ## Check the nr of axts
  checkIdentical(length(axtHg19DanRer7), 133L)

  ## Check the subAxt on target
  foo <- subAxt(axtHg19DanRer7, chr="chr11", start=31500113, end=31500120, 
                select="target")
  checkIdentical(as.character(targetSeqs(foo)), "AAATGCAG")
  checkIdentical(as.character(querySeqs(foo)), "GAGTGC-T")


}

