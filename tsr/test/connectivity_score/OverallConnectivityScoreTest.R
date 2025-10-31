library(testthat)
source("modules/config.R")
source("modules/connectivity_score/OverallConnectivityScore.R")
source("modules/connectivity_score/ConnectivityScoreByGeneSelectionBinChen.R")

# setup
sut <- OverallConnectivityScore$new(ConnectivityScoreByGeneSelectionBinChen$new())

# given
cs_parameters <- list()
cs_parameters$disease_name <- "ESCA"
cs_parameters$disease_most_significant_genes_filename <- "test/connectivity_score/data/ESCA_signature.Rds"
cs_parameters$drugs_filter <- NA
cs_parameters$perturbation_time_list <- c("6h", "24h", "6h_24h")
cs_parameters$gene_selection_strategy <- disease_signature_cfg$gene_select_strategy_150_most_significant
cs_parameters$n_permutations <- 10
expected <- readRDS("test/connectivity_score/data/overall_connectivity_score_expected.Rds")

# when
result <- sut$compute(cs_parameters)
result[, 6] <- NULL
#saveRDS(result, "test/connectivity_score/data/overall_connectivity_score_expected.Rds")
#write.table(result, "test/connectivity_score/data/overall_connectivity_score_expected.csv", sep = ";", row.names = F)

# then
test_that("overall_connectivity_score_test", {
  expect_identical(result, expected)
}
)