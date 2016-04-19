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
  data(axtHg19DanRer7)
  
  ## Test selection on target sequence
  axt <- subAxt(axtHg19DanRer7, chr="chr11", start=31500000L, end=32500000L,
                select="target")
  expect_identical(length(axt), 94L)
  
  ## Test selection on query sequence
  searchGRanges <- GRanges(seqnames="chr25",
                           ranges=IRanges(start=15559655,
                                          end=15575192),
                           strand="+")
  axt <- subAxt(axtHg19DanRer7, searchGRanges, select="query", 
                qSize=c("chr25"=38499472L))
  expect_identical(length(axt), 5L)
})

