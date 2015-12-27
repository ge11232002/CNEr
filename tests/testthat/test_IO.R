test_that("test_axtInfo", {
  axtFilesHg19DanRer7 <- file.path(system.file("extdata", package="CNEr"), 
                                   "hg19.danRer7.net.axt")
  ans <- axtInfo(axtFilesHg19DanRer7)
  expect_identical(133L, length(ans))
  expect_identical(c(479L, 316L, 136L, 103L, 92L), ans[1:5])
  }
)



