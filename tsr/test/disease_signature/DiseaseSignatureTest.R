library(testthat)

source("modules/disease_signature/DiseaseSignature.R")

# setup
sut <- DiseaseSignature$new()

# given
gene_experiments_data <- readRDS("test/disease_signature/data/gene_experiments_data.Rds")
expected <- readRDS("test/disease_signature/data/expected.Rds")

# when
result <- sut$compute(gene_experiments_data, F)
result$p.value <- NULL
result$adj.p.value <- NULL
#saveRDS(result, "test/disease_signature/data/expected.Rds")

# then
test_that("DiseaseSignatureTest", {
  expect_identical(result, expected)
}
)