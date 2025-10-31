library(testthat)
source("modules/config.R")
source("modules/drug_signature/OverallMetanalysis6h24hDrugSignature.R")

# setup
sut <- OverallMetanalysis6h24hDrugSignature$new("test/drug_signature/data/signatures/genes/", "test/drug_signature/data/signatures/metanalysis/")

# given
gene_list_filename <- "test/drug_signature/data/gene_list.csv"
expected_A2M <- readRDS("test/drug_signature/data/signatures/metanalysis/A2M_2_metanalysis-expected.Rds")
expected_EPHA3 <- readRDS("test/drug_signature/data/signatures/metanalysis/EPHA3_2042_metanalysis-expected.Rds")

# when
no_output <- sut$compute(gene_list_filename)
result_A2M <- readRDS("test/drug_signature/data/signatures/metanalysis/A2M_2_metanalysis.Rds")
result_EPHA3 <- readRDS("test/drug_signature/data/signatures/metanalysis/EPHA3_2042_metanalysis.Rds")

# then
test_that("OverallMetanalysis6h24hDrugSignatureTest", {
  expect_equal(result_A2M, expected_A2M)
  expect_equal(result_EPHA3, expected_EPHA3)
}
)