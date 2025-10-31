library(testthat)

source("modules/disease_signature/DiseaseSignature.R")

# setup
sut <- DiseaseSignature$new()

# given
gene_experiments_data <- readRDS("test/disease_signature/data/disease_signature_filter_by_protein_coding_input.Rds")
expected <- readRDS("test/disease_signature/data/disease_signature_filter_by_protein_coding_expected.Rds")

# when
result <- sut$compute(gene_experiments_data)
result$p.value <- NULL
result$adj.p.value <- NULL
saveRDS(result, "test/disease_signature/data/disease_signature_filter_by_proteing_coding_expected.Rds")

# then
test_that("DiseaseSignatureFilterByProteinCodingTest", {
  expect_identical(result, expected)
}
)