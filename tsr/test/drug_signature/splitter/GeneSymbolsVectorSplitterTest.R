library(testthat)
source("modules/drug_signature/splitter/GeneSymbolsVectorSplitter.R")

# setup
sut <- GeneSymbolsVectorSplitter$new()

# given
filename <- "output/disease_signature/als_1/als_1_150_most_significant_genes.csv"

# when
result <- sut$split(filename, 7)

#then
test_that("GeneSymbolsVectorSplitterTest", {
  for (i in 1:6) {
    expect_equal(length(result[[i]]), 22)
  }
  expect_equal(length(result[[7]]), 18)
}
)