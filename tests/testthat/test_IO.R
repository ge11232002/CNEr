test_that("test_axtInfo", {
  bedHg19Fn <- file.path(system.file("extdata", package="CNEr"), 
                         "filter_regions.hg19.bed")
  bedHg19 <- readBed(bedHg19Fn)
  ## Chekc the number of items
  expect_identical(length(bedHg19), 5574L)
  ## Check the first 5 ranges
  expect_identical(ranges(bedHg19[1:5]),
                   IRanges(start=c(30000001, 30000126, 30001660, 
                                   30003969, 30004481),
                           end=c(30000106, 30000274, 30003939, 30004060,
                                 30004844)))
  }
)

test_that("test_axtInfo", {
  axtFilesHg19DanRer7 <- file.path(system.file("extdata", package="CNEr"), 
                                   "hg19.danRer7.net.axt")
  axtInfoHg19DanRer7 <- axtInfo(axtFilesHg19DanRer7)
  ## Check the number of alignments
  expect_identical(length(axtInfoHg19DanRer7), 133L)
  ## Check the first five widths of alignments
  expect_identical(axtInfoHg19DanRer7[1:5], c(479L, 316L, 136L, 103L, 92L))
  }
)

test_that("test_readAxt", {
  axtFilesHg19DanRer7 <- file.path(system.file("extdata", package="CNEr"), 
                                   "hg19.danRer7.net.axt")
  axtHg19DanRer7 <- readAxt(axtFilesHg19DanRer7)
  ## Check the number of alignments
  expect_identical(length(axtHg19DanRer7), 133L)
  ## Check the target ranges
  library(GenomicRanges)
  ansTargetsRanges <- GRanges(seqnames=c("chr11"),
                              ranges=IRanges(start=c(31082021, 31082562, 
                                                     31085908,
                                                     31087045, 31090791),
                                             end=c(31082458, 31082862, 
                                                   31086035,
                                                   31087145, 31090882)),
                              strand="+")
  expect_identical(seqnames(axtHg19DanRer7@targetRanges[1:5]),
                   seqnames(ansTargetsRanges))
  expect_identical(ranges(axtHg19DanRer7@targetRanges[1:5]),
                   ranges(ansTargetsRanges))
  expect_identical(strand(axtHg19DanRer7@targetRanges[1:5]),
                   strand(ansTargetsRanges))
  ## Check the query ranges
  ansQueryRanges <- GRanges(seqnames=c("chr25", "chr25", "chr7", 
                                       "chr25", "chr7"),
                            ranges=IRanges(start=c(15030563, 15036393, 
                                                   60364309, 15037012, 
                                                   60364679),
                                           end=c(15031009, 15036688, 60364442,
                                                 15037110, 60364769)),
                            strand=c("+", "+", "-", "+", "-")) 
  expect_identical(as.character(seqnames(axtHg19DanRer7@queryRanges[1:5])),
                   as.character(seqnames(ansQueryRanges)))
  expect_identical(ranges(axtHg19DanRer7@queryRanges[1:5]),
                   ranges(ansQueryRanges))
  expect_identical(strand(axtHg19DanRer7@queryRanges[1:5]), 
                   strand(ansQueryRanges))
  
  ## Check the score
  expect_identical(axtHg19DanRer7@score[1:5],
                   c(246L, 4422L, 5679L, 1743L, 1556L))
  
  library(Biostrings)
  ## Check the target seqs
  expect_equivalent(subseq(axtHg19DanRer7@targetSeqs[1:5], start=1, end=5),
                    DNAStringSet(c("ATTTT", "GGGAA", "GGGCT", "TTGTG", "TATTC"))
                    )
  ## Check the query seqs
  expect_equivalent(subseq(axtHg19DanRer7@querySeqs[1:5], start=1, end=5),
                    DNAStringSet(c("ATTTA", "GGAAA", "GGGCT", "TTAAA", "TATTT"))
                    )
  }
)

