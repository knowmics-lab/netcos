library(testthat)
source("modules/config.R")
source("modules/connectivity_score/ConnectivityScoreByGeneSelectionBinChen.R")
source("modules/connectivity_score/ConnectivityScoreByPerturbationTime.R")

# setup
sut <- ConnectivityScoreByPerturbationTime$new(ConnectivityScoreByGeneSelectionBinChen$new())

# given
#drug_by_gene_signatures <- read.table("test/connectivity_score/data/overall_meta_analysis_drug_by_gene_signature_6h_24h.csv", sep = ";", header = T)
#drug_by_gene_signatures$t.value_6h<-as.numeric(drug_by_gene_signatures$t.value_6h)
#saveRDS(drug_by_gene_signatures,"test/connectivity_score/data/overall_meta_analysis_drug_by_gene_signature_6h_24h.Rds")

disease_name <- "BLCA"
disease_signature_filename <- "test/connectivity_score/data/disease_signature_TCGA.Rds"
drug_signatures_filename <- "test/connectivity_score/data/overall_meta_analysis_drug_by_gene_signature_6h_24h.Rds"
drugs_filter <- NA
perturbation_time <- "6h"
expected <- readRDS("test/connectivity_score/data/connectivity_score_by_perturbation_time_expected.Rds")

# when
result <- sut$compute(disease_name, disease_signature_filename, drug_signatures_filename, drugs_filter, perturbation_time, disease_signature_cfg$gene_select_strategy_150_most_significant, n_permutations = 10)
result[, 6] <- NULL
#saveRDS(result, "test/connectivity_score/data/connectivity_score_by_perturbation_time_expected.Rds")
#write.table(result, "test/connectivity_score/data/connectivity_score_by_perturbation_time_expected.csv", sep = ";", row.names = F)

# then
test_that("ConnectivityScoreByStrategyTest", {
  expect_identical(result, expected)
}
)