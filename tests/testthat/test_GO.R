test_that("test_addAncestorGO", {
  library(GO.db)
  go <- list(c("GO:0005215", "GO:0006810", "GO:0016020"), "GO:0016579")
  newGO <- addAncestorGO(go)
  expect_identical(lengths(newGO), c(8, 17))
})