library(testthat)
source("modules/connectivity_score/ConnectivityScoreByGeneSelectionBinChen.R")
source("modules/connectivity_score/mapper/DiseaseSignatureMapper.R")

# setup
sut <- ConnectivityScoreByGeneSelectionBinChen$new()
diseaseSignatureMapper <- DiseaseSignatureMapper$new()

# given
drug_signatures <- readRDS("test/connectivity_score/data/6h_drug_signatures_LINCS.Rds")
colnames(drug_signatures)[5] <- "t.value"
disease_signature <- readRDS("test/connectivity_score/data/BLCA_signature.Rds")
disease_signature <- disease_signature[order(abs(disease_signature$t.value), decreasing = T),]
disease_signature <- disease_signature[1:50,]
disease_signature <- diseaseSignatureMapper$map(disease_signature)
expected <- readRDS("test/connectivity_score/data/connectivity_score_by_disease_BLCA_expected.Rds")

# when
result <- sut$compute("BLCA", disease_signature, drug_signatures, "50_most_significant", "6h", 10)
result$p.value <- NULL
#saveRDS(result, "test/connectivity_score/data/connectivity_score_by_disease_BLCA_expected.Rds")
#write.table(result, "test/connectivity_score/data/connectivity_score_by_disease_BLCA_expected.csv", sep = ";", row.names = F)

# then
test_that("ConnectivityScoreByGeneSelectionBLCATest", {
  expect_identical(result, expected)
}
)

