library(testthat)
source("modules/connectivity_score/loader/DrugSignaturesLoader.R")

# setup
sut <- DrugSignaturesLoader$new()

# given
drug_signatures_filename <- "test/connectivity_score/data/drug_by_gene_signature_6h_24h_metanalysis.Rds"
perturbation_time <- "6h_24h"
disease_genes <- c("DDR1", "FAU")
expected <- readRDS("test/connectivity_score/data/load_and_setup_drug_signatures_6h_24h_expected.Rds")

# when
result <- sut$load(drug_signatures_filename, disease_genes, perturbation_time)
#saveRDS(result, "test/connectivity_score/data/load_and_setup_drug_signatures_6h_24h_expected.Rds")
#write.table(result, "test/connectivity_score/data/load_and_setup_drug_signatures_6h_24h_expected.csv", sep = ";", row.names = F)


# then
test_that("LoadDrugSignatures6h24hTest", {
  expect_identical(result, expected)
}
)