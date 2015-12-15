test_that("test_reverseCigar", {
    cigar <- "16I20M17D"
    expect_identical("17D20M16I", reverseCigar(cigar))
    cigar <- "17D20M16I"
    expect_identical("16I20M17D", reverseCigar(cigar))
    cigar <- "20M17D16I"
    expect_identical("16I17D20M", reverseCigar(cigar))
})

