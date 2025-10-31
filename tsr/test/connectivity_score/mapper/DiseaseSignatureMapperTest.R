library(testthat)
source("modules/connectivity_score/mapper/DiseaseSignatureMapper.R")
source("modules/connectivity_score/loader/DiseaseSignatureLoader.R")
source("modules/config.R")

# setup
sut <- DiseaseSignatureMapper$new()
diseaseSignatureLoader <- DiseaseSignatureLoader$new()

# given
disease_signature_filename <- gsub(".csv", ".Rds", ipf_disease_signature_cfg$ipf_150_most_significant_genes_filename)
disease_signature <- diseaseSignatureLoader$load(disease_signature_filename)
disease_signature <- disease_signature[order(disease_signature$gene),]
expected <- readRDS("test/connectivity_score/mapper/data/DiseaseSignatureMapper.Rds")

# when
result <- result <- sut$map(disease_signature)
#write.table(result, "test/connectivity_score/mapper/data/DiseaseSignatureMapper.csv", sep = ";", dec = ",", row.names = T)
#saveRDS(result, "test/connectivity_score/mapper/data/DiseaseSignatureMapper.Rds")

# then
test_that("DiseaseSignatureMapperTest", {
  expect_identical(result, expected)
}
)