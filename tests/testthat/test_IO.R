test_that("test_readBed", {
  bedHg19Fn <- file.path(system.file("extdata", package="CNEr"), 
                         "filter_regions.hg38.bed")
  bedHg19 <- readBed(bedHg19Fn)
  ## Chekc the number of items
  expect_identical(length(bedHg19), 413L)
  ## Check the first 5 ranges
  expect_identical(start(bedHg19[1:5]),
                   as.integer(c(26817474, 156772018, 156772160, 
                                156772495, 156774135)))
  }
)

test_that("test_axtInfo", {
  axtFilesHg19DanRer7 <- file.path(system.file("extdata", package="CNEr"), 
                                   "hg38.danRer10.net.axt")
  axtInfoHg19DanRer7 <- axtInfo(axtFilesHg19DanRer7)
  ## Check the number of alignments
  expect_identical(length(axtInfoHg19DanRer7), 50L)
  ## Check the first five widths of alignments
  expect_identical(axtInfoHg19DanRer7[1:5], c(359L, 398L, 961L, 673L, 55L))
  }
)

test_that("test_readAxt", {
  axtFilesHg19DanRer7 <- file.path(system.file("extdata", package="CNEr"), 
                                   "hg38.danRer10.net.axt")
  axtHg19DanRer7 <- readAxt(axtFilesHg19DanRer7)
  ## Check the number of alignments
  expect_identical(length(axtHg19DanRer7), 50L)
  ## Check the target ranges
  library(GenomicRanges)
  ansTargetsRanges <- GRanges(seqnames=c("chr1", "chr1", "chr10", "chr13", 
                                         "chr13"),
                              ranges=IRanges(start=c(148165963, 222131480, 
                                                     65322021,
                                                     94750629, 94966940),
                                             end=c(148166304, 222131835, 
                                                   65322919,
                                                   94751259, 94966991)),
                              strand="+")
  expect_identical(as.character(seqnames(targetRanges(axtHg19DanRer7)[1:5])),
                   as.character(seqnames(ansTargetsRanges)))
  expect_identical(ranges(targetRanges(axtHg19DanRer7)[1:5]),
                   ranges(ansTargetsRanges))
  expect_identical(strand(targetRanges(axtHg19DanRer7)[1:5]),
                   strand(ansTargetsRanges))
  ## Check the query ranges
  ansQueryRanges <- GRanges(seqnames=c("chr6"),
                            ranges=IRanges(start=c(25825774, 24819722, 
                                                   26227619, 26600600, 
                                                   26745445),
                                           end=c(25826099, 24820074, 26228468,
                                                 26601208, 26745499)),
                            strand=c("+")) 
  expect_identical(as.character(seqnames(queryRanges(axtHg19DanRer7)[1:5])),
                   as.character(seqnames(ansQueryRanges)))
  expect_identical(ranges(queryRanges(axtHg19DanRer7)[1:5]),
                   ranges(ansQueryRanges))
  expect_identical(strand(queryRanges(axtHg19DanRer7)[1:5]), 
                   strand(ansQueryRanges))
  
  ## Check the score
  expect_identical(score(axtHg19DanRer7)[1:5],
                   c(5310L, 221L, 946L, 3302L, 711L))
  
  library(Biostrings)
  ## Check the target seqs
  expect_equivalent(subseq(targetSeqs(axtHg19DanRer7)[1:5], start=1, end=5),
                    DNAStringSet(c("TCTCT", "CTTAA", "TTGCA", "GCCTT", 
                                   "ATTCA"))
                    )
  ## Check the query seqs
  expect_equivalent(subseq(querySeqs(axtHg19DanRer7)[1:5], start=1, end=5),
                    DNAStringSet(c("TTTCT", "TTCGA", "TTGTT", "GCTTT", 
                                   "AGCCA"))
                    )
  }
)

test_that("test_writeAxt", {
  axtFilesHg19DanRer7 <- file.path(system.file("extdata", package="CNEr"), 
                                   "hg38.danRer10.net.axt")
  axtHg19DanRer7 <- readAxt(axtFilesHg19DanRer7)
  
  ## Check we can output the axt
  expect_silent(writeAxt(axtHg19DanRer7, con=tempfile()))
  
})

test_that("test_read.rmMask.GRanges", {
  fn <- system.file("extdata", "ce2chrM.fa.out", package="IRanges")
  ans <- read.rmMask.GRanges(fn)
  
  ## Check the number of repeats
  expect_identical(length(ans), 28L)
})

test_that("test_readCNERangesFromSQLite", {
  dbName <- file.path(system.file("extdata", package="CNEr"),
                      "danRer10CNE.sqlite")
  cneGRangePairs <- readCNERangesFromSQLite(dbName=dbName, 
                                            tableName="danRer10_hg38_45_50")
  expect_equal(length(cneGRangePairs), 3209L)
})

test_that("test_read.rmskFasta", {
  library(GenomicRanges)
  fn <- file.path(system.file("extdata", package="CNEr"),
                  "rmsk.fa")
  ans <- read.rmskFasta(fn)
  expect_equal(ans, GRanges(seqnames=c("chr1", "chr2"),
                            ranges=IRanges(start=c(6,3),
                                           end=c(7,5)),
                            strand="+"))
})