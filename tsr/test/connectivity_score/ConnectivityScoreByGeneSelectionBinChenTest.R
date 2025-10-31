library(testthat)
source("modules/connectivity_score/ConnectivityScoreByGeneSelectionBinChen.R")
source("modules/connectivity_score/mapper/DiseaseSignatureMapper.R")

# setup
sut <- ConnectivityScoreByGeneSelectionBinChen$new()

# given
drug_by_gene_signatures <- readRDS("test/connectivity_score/data/drug_by_gene_signatures_6h.Rds")
disease_signature <- data.frame(
  estimate = c(-1.101117, 0.91904, 0.202952)
)
rownames(disease_signature) <- c("DDR1", "PAX8", "FAU")
expected <- readRDS("test/connectivity_score/data/connectivity_score_by_disease.Rds")

# when
result <- sut$compute("als", disease_signature, drug_by_gene_signatures, "50_most_significant", "6h", 10)
result$p.value <- NULL
#saveRDS(result, "test/connectivity_score/data/connectivity_score_by_disease.Rds")
#write.table(result, "test/connectivity_score/data/connectivity_score_by_disease.csv", sep = ";", row.names = F)

# then
test_that("ConnectivityScoreByGeneSelectionTest", {
  expect_identical(result, expected)
}
)
