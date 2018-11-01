test_that("test_ceScan", {
  axtFnHg38DanRer10 <- file.path(system.file("extdata", package="CNEr"),
                                 "hg38.danRer10.net.axt")
  axtHg38DanRer10 <- readAxt(axtFnHg38DanRer10)
  axtFnDanRer10Hg38 <- file.path(system.file("extdata", package="CNEr"),
                                 "danRer10.hg38.net.axt")
  axtDanRer10Hg38 <- readAxt(axtFnDanRer10Hg38)
  bedHg38Fn <- file.path(system.file("extdata", package="CNEr"),
                         "filter_regions.hg38.bed")
  bedHg38 <- readBed(bedHg38Fn)
  bedDanRer10Fn <- file.path(system.file("extdata", package="CNEr"),
                             "filter_regions.danRer10.bed")
  bedDanRer10 <- readBed(bedDanRer10Fn)
  library(BSgenome.Drerio.UCSC.danRer10)
  library(BSgenome.Hsapiens.UCSC.hg38)
  qSizesHg38 <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
  qSizesDanRer10 <- seqinfo(BSgenome.Drerio.UCSC.danRer10)
  
  windows <- c(50L, 50L, 50L)
  identities <- c(45L, 48L, 49L)
  CNEHg38DanRer10 <- ceScan(x=axtHg38DanRer10, tFilter=bedHg38,
                            qFilter=bedDanRer10,
                            tSizes=qSizesHg38, qSizes=qSizesDanRer10,
                            window=windows, identity=identities)
  # check the names of CNE list
  expect_identical(names(CNEHg38DanRer10), c("45_50", "48_50", "49_50"))
  
  # check the CNEs in "45_50" set
  expect_identical(length(CNEHg38DanRer10[["45_50"]]), 4L)
  
  targetCNEs <- GRanges(seqnames=c("chr3"),
                        ranges=IRanges(start=c(137523038, 137523122,
                                               137269941, 137264717),
                                       end=c(137523102, 137523187,
                                             137269991, 137265124)),
                        strand="+", seqinfo=qSizesHg38)
  expect_identical(first(CNEHg38DanRer10[["45_50"]]), targetCNEs)
  
  queryCNEs <- GRanges(seqnames=c("chr6"),
                       ranges=IRanges(start=c(26638744, 26638826,
                                              26746284, 26745047),
                                      end=c(26638808, 26638891,
                                            26746334, 26745455)),
                       strand="+", seqinfo=qSizesDanRer10)
  expect_identical(last(CNEHg38DanRer10[["45_50"]]), queryCNEs)
  
  library(S4Vectors)
  mcolsCNEs <- DataFrame(score=c(87.69, 89.39, 90.20, 95.11),
                        cigar=c("65M", "66M", "51M", "218M1D190M"))
  expect_identical(mcols(CNEHg38DanRer10[["45_50"]]), mcolsCNEs)
}
)