test_that("test_CNE", {
  library(GenomicRanges)
  CNE12 <- GRangePairs(first=GRanges(seqnames=c("chr13", "chr4", "chr4"),
                                     ranges=IRanges(start=c(71727138,150679343,
                                                            146653164),
                                                    end=c(71727224, 150679400,
                                                          146653221)),
                                     strand="+"),
                       last=GRanges(seqnames=c("chr1"),
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
                       last=GRanges(seqnames=c("chr13", "chr4", "chr4"),
                                    ranges=IRanges(start=c(71727138,150679343,
                                                           146653164),
                                                   end=c(71727224, 150679400,
                                                         146653221)),
                                    strand="+")
  )
  
  # Test the validity
  ## test window < identity
  expect_error(CNE(assembly1Fn="hg38.2bit", assembly2Fn="danRer10.2bit", 
                   window=49L,
                   identity=50L, CNE12=CNE12, CNE21=CNE21,
                   CNEMerged=CNE12, CNEFinal=CNE12, aligner="blat",
                   smoothingWindow1=300L, smoothingWindow2=300L))
  ## test number of assembly, aligner >= 2
  expect_error(CNE(assembly1Fn=c("hg38.2bit", "galGal4.2bit"), 
                   assembly2Fn="danRer10.2bit", window=50L,
                   identity=50L, CNE12=CNE12, CNE21=CNE21,
                   CNEMerged=CNE12, CNEFinal=CNE12, aligner="blat",
                   smoothingWindow1=300L, smoothingWindow2=300L))
  expect_error(CNE(assembly1Fn="hg38.2bit",
                   assembly2Fn=c("danRer10.2bit", "tetNig2.bit"), window=50L,
                   identity=50L, CNE12=CNE12, CNE21=CNE21,
                   CNEMerged=CNE12, CNEFinal=CNE12, aligner="blat",
                   smoothingWindow1=300L, smoothingWindow2=300L))
  expect_error(CNE(assembly1Fn="hg38.2bit", assembly2Fn="danRer10.2bit",
                   window=50L,
                   identity=50L, CNE12=CNE12, CNE21=CNE21,
                   CNEMerged=CNE12, CNEFinal=CNE12, aligner=c("blat", "bwa"),
                   smoothingWindow1=300L, smoothingWindow2=300L))
  ## test smoothing window
  expect_error(CNE(assembly1Fn="hg38.2bit", assembly2Fn="danRer10.2bit",
                   window=50L,
                   identity=50L, CNE12=CNE12, CNE21=CNE21,
                   CNEMerged=CNE12, CNEFinal=CNE12, aligner=c("blat"),
                   smoothingWindow1=5L, smoothingWindow2=300L))
  expect_error(CNE(assembly1Fn="hg38.2bit", assembly2Fn="danRer10.2bit",
                   window=50L,
                   identity=50L, CNE12=CNE12, CNE21=CNE21,
                   CNEMerged=CNE12, CNEFinal=CNE12, aligner=c("blat"),
                   smoothingWindow1=2000L, smoothingWindow2=300L))
  
  # Test the getter
  cne <- CNE(assembly1Fn="hg38.2bit", assembly2Fn="danRer10.2bit",
             window=50L, identity=50L,
             CNE12=CNE12, CNE21=CNE21, CNEMerged=CNE12, CNEFinal=CNE12,
             aligner="blat", smoothingWindow1=300L, smoothingWindow2=300L)
  #expect_identical(assembly1(cne), "hg38.2bit")
  #expect_identical(assembly2(cne), "danRer10.2bit")
  expect_identical(CNE12(cne), CNE12)
  expect_identical(CNE21(cne), CNE21)
  expect_identical(thresholds(cne), "50_50")
  expect_identical(CNEMerged(cne), CNE12)
  expect_identical(CNEFinal(cne), CNE12)
}
)