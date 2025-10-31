library(testthat)
source("modules/drug_signature/splitter/GeneSymbolsVectorSplitter.R")
source("modules/drug_signature/builder/DrugSignatureInputChunksBuilder.R")

# setup
geneSymbolsVectorSplitter <- GeneSymbolsVectorSplitter$new()
sut <- DrugSignatureInputChunksBuilder$new()

# given
filename <- "output/disease_signature/als_1/als_1_150_most_significant_genes.csv"
gene_symbols <- geneSymbolsVectorSplitter$split(filename, 4)
BLAS_num_threads_list <- c(1, 4, 3, 5, 9, 11, 7, 8)

# when
result <- sut$build(gene_symbols, c("6", "24"), BLAS_num_threads_list)

#then
test_that("DrugSignatureSplitInputTest", {
  for (i in 1:8) {
    expect_equal(result[[i]]$BLAS_num_threads, BLAS_num_threads_list[i])
  }
  for (i in 1:3) {
    expect_equal(length(result[[i]]$gene_symbols), 38)
    expect_equal(length(result[[i + 4]]$gene_symbols), 38)
  }
  expect_equal(length(result[[4]]$gene_symbols), 36)
  expect_equal(length(result[[8]]$gene_symbols), 36)
  for (i in 1:4) {
    expect_equal(result[[i]]$perturbation_time, "6")
    expect_equal(result[[i + 4]]$perturbation_time, "24")
  }
}
)