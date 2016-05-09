test_that("test_makeGRBs", {
  library(TxDb.Drerio.UCSC.danRer10.refGene)
  refGenesDanRer10 <- genes(TxDb.Drerio.UCSC.danRer10.refGene)
  ancoraCNEsFns <- file.path(system.file("extdata", package="CNEr"), 
                             c("cne2wBf_cypCar1_danRer10_100_100",
                               "cne2wBf_cteIde1_danRer10_100_100",
                               "cne2wBf_AstMex102_danRer10_48_50"))
  cneList <- do.call(GRangesList, 
                     lapply(ancoraCNEsFns, readAncora, assembly="danRer10"))
  names(cneList) <- c("Common carp", "Grass carp", "Blind cave fish")
  seqlengths(cneList) <- seqlengths(TxDb.Drerio.UCSC.danRer10.refGene)[
    names(seqlengths(cneList))]
  danRer10GRBs <- makeGRBs(cneList, winSize=200, genes=refGenesDanRer10,
                           ratio=1.2)
  expect_identical(length(danRer10GRBs), 502L)
}
)