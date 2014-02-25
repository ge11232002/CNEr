

test_reverseCigar <- function(){
  cigar <- "16I20M17D"
  checkIdentical(reverseCigar(cigar), "17D20M16I")
  cigar <- "17D20M16I" 
  checkIdentical(reverseCigar(cigar), "16I20M17D")
  cigar <- "20M17D16I"
  checkIdentical(reverseCigar(cigar), "16I17D20M")

}

