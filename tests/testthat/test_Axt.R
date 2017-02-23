test_that("test_Axt", {
  library(GenomicRanges)
  library(Biostrings)
  targetRanges <- GRanges(seqnames=c("chr1", "chr1", "chr2", "chr3"),
                          ranges=IRanges(start=c(1, 20, 2, 3),
                                         end=c(10, 25, 10, 10)),
                          strand="+")
  targetSeqs <- DNAStringSet(c("ATTTTATGTG", "GGGAAG", "GGGCTTTTG",
                               "TTGTGTAG"))
  queryRanges <- GRanges(seqnames=c("chr1", "chr10", "chr10", "chr20"),
                         ranges=IRanges(start=c(1, 25, 50, 5),
                                        end=c(10, 30, 58, 12)),
                         strand="+")
  querySeqs <- DNAStringSet(c("ATTTAAAGTG", "GGAAAA", "GGGCTCTGG", "TTAAATAA"))
  score <- c(246L, 4422L, 5679L, 1743L)
  symCount <- c(10L, 6L, 9L, 8L)
  axt <- Axt(targetRanges=targetRanges, targetSeqs=targetSeqs, 
             queryRanges=queryRanges, querySeqs=querySeqs, 
             score=score, symCount=symCount)
  
  # Test the validity
  ## test the different lengths
  expect_error(Axt(targetRanges=targetRanges, targetSeqs=targetSeqs, 
                   queryRanges=queryRanges, querySeqs=querySeqs, 
                   score=score, symCount=symCount[1:3]))
  expect_error(Axt(targetRanges=targetRanges, targetSeqs=targetSeqs, 
                   queryRanges=queryRanges, querySeqs=querySeqs, 
                   score=score[1:3], symCount=symCount))
  expect_error(Axt(targetRanges=targetRanges[1:3], targetSeqs=targetSeqs[1:3], 
                   queryRanges=queryRanges, querySeqs=querySeqs, 
                   score=score, symCount=symCount))
  expect_error(Axt(targetRanges=targetRanges, targetSeqs=targetSeqs,
                   queryRanges=queryRanges[1:3], querySeqs=querySeqs[1:3],
                   score=score, symCount=symCount))
  ## test the inconsistent widths
  expect_error(Axt(targetRanges=targetRanges, targetSeqs=targetSeqs,
                   queryRanges=queryRanges, querySeqs=querySeqs,
                   score=score, symCount=c(10L, 6L, 9L, 9L)))
  expect_error(Axt(targetRanges=targetRanges,
                   targetSeqs= DNAStringSet(c("ATTTTATGTG", "GGGAAG",
                                              "GGGCTTTTG","TTGTGTA")),
                   queryRanges=queryRanges, querySeqs=querySeqs,
                   score=score, symCount=symCount))
  
  # Test the getters and setters
  ## test length
  expect_identical(length(axt), 4L)
  ## test targetRanges, targetSeqs
  expect_equivalent(targetRanges(axt), targetRanges)
  expect_identical(targetSeqs(axt), targetSeqs)
  ## test queryRanges, querySeqs
  expect_equivalent(queryRanges(axt), queryRanges)
  expect_identical(querySeqs(axt), querySeqs)
  ## test score, symCount
  expect_identical(score(axt), score)
  expect_identical(symCount(axt), symCount)
  
  # TODO: Test matchDistr
  ## matchDistr
  
  # Test Vector methods
  expect_identical(length(axt[1:2]), 2L)
  
  # test combing
  expect_identical(length(c(axt, axt)), 8L)
}
)

test_that("test_subAxt", {
  library(GenomicRanges)
  library(rtracklayer)
  tAssemblyFn <- file.path(system.file("extdata",
                                       package="BSgenome.Hsapiens.UCSC.hg38"),
                           "single_sequences.2bit")
  qAssemblyFn <- file.path(system.file("extdata",
                                       package="BSgenome.Drerio.UCSC.danRer10"),
                           "single_sequences.2bit")
  axtFilesHg38DanRer10 <- file.path(system.file("extdata", package="CNEr"),
                                    "hg38.danRer10.net.axt")
  axtHg38DanRer10 <- readAxt(axtFilesHg38DanRer10, tAssemblyFn, qAssemblyFn)

  ## Test selection on target sequence
  axt <- subAxt(axtHg38DanRer10, chr="chr1", start=148165963L, end=222131835L,
                select="target")
  expect_identical(length(axt), 2L)
  
  ## Test selection on query sequence
  searchGRanges <- GRanges(seqnames="chr6",
                           ranges=IRanges(start=25825774,
                                          end=26745499),
                           strand="+")
  axt <- subAxt(axtHg38DanRer10, searchGRanges, select="query")
  expect_identical(length(axt), 9L)
})

# test_that("test_summaryAxt", {
#   axtFilesHg38DanRer10 <- file.path(system.file("extdata", package="CNEr"),
#                                     "hg38.danRer10.net.axt")
#   axtHg38DanRer10 <- readAxt(axtFilesHg38DanRer10)
#   ans <- summary(axtHg38DanRer10)
# #  expect_equal(as.integer(ans),
# #               c(907L, 5012L, 1133L, 3274L, 1801L, 1606L, 3078L))
# }
# )

test_that("test_fixCoordinates", {
  axtFnDanRer10Hg38 <- file.path(system.file("extdata", package="CNEr"),
                                 "danRer10.hg38.net.axt")
  qAssemblyFn <- file.path(system.file("extdata",
                                       package="BSgenome.Hsapiens.UCSC.hg38"),
                           "single_sequences.2bit")
  tAssemblyFn <- file.path(system.file("extdata",
                                       package="BSgenome.Drerio.UCSC.danRer10"),
                           "single_sequences.2bit")
  axtDanRer10Hg38 <- readAxt(axtFnDanRer10Hg38, tAssemblyFn=tAssemblyFn,
                             qAssemblyFn=qAssemblyFn)
  fixedAxt <- fixCoordinates(axtDanRer10Hg38)
  # Test the restore
  expect_identical(fixCoordinates(fixedAxt), axtDanRer10Hg38)
  # Test that targetRanges are intact
  expect_identical(targetRanges(fixedAxt), targetRanges(axtDanRer10Hg38))
  # Test the fixed queryRanges coordinates.
  queryRangesTest <- queryRanges(fixedAxt)[1:3]
  expect_identical(ranges(queryRangesTest), IRanges(start=c(12578221, 121302901,
                                                            12577855),
                                                    end=c(12578959, 121303067,
                                                          12578037))
                   )
}
)

test_that("test_subAxt", {
  tAssemblyFn <- file.path(system.file("extdata",
                                       package="BSgenome.Drerio.UCSC.danRer10"),
                           "single_sequences.2bit")
  qAssemblyFn <- file.path(system.file("extdata",
                                       package="BSgenome.Hsapiens.UCSC.hg38"),
                           "single_sequences.2bit")
  axtFn <- file.path(system.file("extdata", package="CNEr"), 
                     "danRer10.hg38.net.axt")
  axt <- readAxt(axtFn, tAssemblyFn, qAssemblyFn)
  
  targetSearch <- GRanges(seqnames=c("chr6"),
                          ranges=IRanges(start=c(24000000, 26900000),
                                         end=c(24060000, 26905000)),
                          strand="+"
  )
  querySearch <- GRanges(seqnames=c("chr7", "chr2"),
                         ranges=IRanges(start=c(12577000, 241262700),
                                        end=c(12579000, 241268600)),
                         strand="+"
  )
  ans <- psubAxt(axt, targetSearch, querySearch)
  expect_identical(ans, axt[c(1,3,347:350)])
}
)

test_that("test_makeAxtTracks", {
  tAssemblyFn <- file.path(system.file("extdata",
                                       package="BSgenome.Drerio.UCSC.danRer10"),
                           "single_sequences.2bit")
  qAssemblyFn <- file.path(system.file("extdata",
                                       package="BSgenome.Hsapiens.UCSC.hg38"),
                           "single_sequences.2bit")
  axtFn <- file.path(system.file("extdata", package="CNEr"), 
                     "danRer10.hg38.net.axt")
  axt <- readAxt(axtFn, tAssemblyFn, qAssemblyFn)
  ans <- makeAxtTracks(axt)
  
  ## Make sure the coordinates are right.
  expect_identical(ans[[1]]$name[1:2], c("chr7:12578221-12578959:-", 
                                         "chr9:121302901-121303067:-"))
})
