library(testthat)
source("modules/connectivity_score/cmap_score_new/RandomConnectivityScoreParallel.R")

# setup
sut <- RandomConnectivityScoreParallel$new()

# given
n_genes_up <- 31
n_genes_down <- 19
n_genes <- 150
n_permutations_1 <- 233
n_permutations_2 <- 10
n_permutations_3 <- 80
n_permutations_4 <- 807

# when
result_1 <- sut$compute(n_genes_up, n_genes_down, n_genes, n_permutations_1)
result_2 <- sut$compute(n_genes_up, n_genes_down, n_genes, n_permutations_2)
result_3 <- sut$compute(n_genes_up, n_genes_down, n_genes, n_permutations_3)
result_4 <- sut$compute(n_genes_up, n_genes_down, n_genes, n_permutations_4)

# then
test_that("ConnectivityScoreRandomDistributionParallelTest", {
  expect_equal(length(result_1), n_permutations_1)
  expect_equal(length(result_2), n_permutations_2)
  expect_equal(length(result_3), n_permutations_3)
  expect_equal(length(result_4), n_permutations_4)
}
)