library(testthat)
source("modules/connectivity_score/GeneSelectionByDiseaseGeneCount.R")

# setup
sut <- GeneSelectionByDiseaseGeneCount$new()

gene_selection <- sut$compute(150)

# when
result1 <- sut$compute(250)
result2 <- sut$compute(101)
result3 <- sut$compute(100)
result4 <- sut$compute(99)
result5 <- sut$compute(19)

# then
test_that("GeneSelectionByDiseaseGeneCountTest", {
  expect_identical(result1, gene_selection)
  expect_identical(result2$count, c(50, 100, 101))
  expect_identical(result3$count, c(50, 100))
  expect_identical(result4$count, c(50, 99))
  expect_identical(result5$count, 19)
}
)