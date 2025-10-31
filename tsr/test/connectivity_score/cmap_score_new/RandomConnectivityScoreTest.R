library(testthat)
source("modules/connectivity_score/cmap_score_new/RandomConnectivityScore.R")

# setup
sut <- RandomConnectivityScore$new()

# given
n_genes_up <- 31
n_genes_down <- 19
n_genes <- 150
n_permutations <- 20

# when
result <- sut$compute(n_genes_up, n_genes_down, n_genes, n_permutations)

# than
test_that("ConnectivityScoreRandomDistributionParallelTest", {
  expect_equal(length(result), n_permutations)
}
)