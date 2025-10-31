library(testthat)
source("modules/drug_signature/Metanalysis6h24hDrugSignatureByGene.R")

# setup
sut <- Metanalysis6h24hDrugSignatureByGene$new("test/drug_signature/data/signatures/genes/")

# given
drug_by_gene_signatures_base_filename <- "A2M_2"
expected <- readRDS("test/drug_signature/data/signatures/metanalysis/A2M_2_metanalysis-expected.Rds")

# when
result <- sut$compute(drug_by_gene_signatures_base_filename, F)
#write.table(result, "test/drug_signature/data/signatures/genes/A2M_2_metanalysis-expected.csv", sep = ";", row.names = F, dec = ",")
#saveRDS(result, "test/drug_signature/data/signatures/genes/A2M_2_metanalysis-expected.Rds")

# then
test_that("Metanalysis6h24hDrugSignatureByGeneTest", {
  expect_equal(result, expected)
}
)