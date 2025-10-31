library(testthat)
source("modules/drug_signature/merger/Metanalysis6h24hMerger.R")

# setup
sut <- Metanalysis6h24hMerger$new()

# given
gene_symbols_filename <- "test/drug_signature/data/metanalysis_gene_list_to_merge.csv"
metanalysis_signature_base_dir <- "test/drug_signature/data/signatures/metanalysis/"
metanalysis_filename <- "test/drug_signature/data/signatures/Metanalysis6h24hMerger.Rds"
expected <- readRDS("test/drug_signature/data/signatures/Metanalysis6h24hMerger-expected.Rds")

# when
no_output <- sut$merge(gene_symbols_filename, metanalysis_signature_base_dir, metanalysis_filename)
result <- readRDS(metanalysis_filename)

# then
test_that("Metanalysis6h24hMergerTest", {
  expect_equal(result, expected)
}
)