test_that("test_CNE", {
  library(GenomicRanges)
  CNE12 <- GRangePairs(first=GRanges(seqnames=c("chr13", "chr4", "chr4"),
                                     ranges=IRanges(start=c(71727138,150679343,
                                                            146653164),
                                                    end=c(71727224, 150679400,
                                                          146653221)),
                                     strand="+"),
                       second=GRanges(seqnames=c("chr1"),
                                      ranges=IRanges(start=c(29854162, 23432387,
                                                             35711077),
                                                     end=c(29854248, 23432444,
                                                           35711134)),
                                      strand="+")
  )
  CNE21 <- GRangePairs(first=GRanges(seqnames=c("chr1"),
                                     ranges=IRanges(start=c(29854162, 23432387,
                                                            35711077),
                                                    end=c(29854248, 23432444,
                                                          35711134)),
                                     strand="+"),
                       second=GRanges(seqnames=c("chr13", "chr4", "chr4"),
                                      ranges=IRanges(start=c(71727138,150679343,
                                                             146653164),
                                                     end=c(71727224, 150679400,
                                                           146653221)),
                                      strand="+")
  )
  assembly1Fn <- file.path(system.file("extdata",
                                    package="BSgenome.Drerio.UCSC.danRer10"),
                        "single_sequences.2bit")
  assembly2Fn <- file.path(system.file("extdata",
                                    package="BSgenome.Hsapiens.UCSC.hg38"),
                        "single_sequences.2bit")
  # Test the validity
  ## test window < identity
  expect_error(CNE(assembly1Fn=assembly1Fn, assembly2Fn=assembly2Fn, 
                   window=49L,
                   identity=50L, CNE12=CNE12, CNE21=CNE21,
                   CNEMerged=CNE12, CNEFinal=CNE12, aligner="blat"))
  ## test number of assembly, aligner >= 2
  expect_error(CNE(assembly1Fn=c(assembly1Fn, assembly2Fn), 
                   assembly2Fn=assembly2Fn, window=50L,
                   identity=50L, CNE12=CNE12, CNE21=CNE21,
                   CNEMerged=CNE12, CNEFinal=CNE12, aligner="blat"))
  expect_error(CNE(assembly1Fn=assembly1Fn,
                   assembly2Fn=c(assembly1Fn, assembly2Fn), window=50L,
                   identity=50L, CNE12=CNE12, CNE21=CNE21,
                   CNEMerged=CNE12, CNEFinal=CNE12, aligner="blat"))
  expect_error(CNE(assembly1Fn=assembly1Fn, assembly2Fn=assembly2Fn,
                   window=50L,
                   identity=50L, CNE12=CNE12, CNE21=CNE21,
                   CNEMerged=CNE12, CNEFinal=CNE12, aligner=c("blat", "bwa")))
  
  # Test the getter
  cne <- CNE(assembly1Fn=assembly1Fn, assembly2Fn=assembly2Fn,
             window=50L, identity=50L,
             CNE12=CNE12, CNE21=CNE21, CNEMerged=CNE12, CNEFinal=CNE12,
             aligner="blat")
  expect_identical(CNE12(cne), CNE12)
  expect_identical(CNE21(cne), CNE21)
  expect_identical(thresholds(cne), "50_50")
  expect_identical(CNEMerged(cne), CNE12)
  expect_identical(CNEFinal(cne), CNE12)
}
)