library(testthat)
source("test/test_config.R")
source("modules/drug_signature/merger/Metanalysis6h24hByGeneSymbolsMerger.R")

# setup
sut <- Metanalysis6h24hByGeneSymbolsMerger$new()

# given
gene_symbols <- c("EPHA3", "A2M")
expected <- readRDS("test/drug_signature/data/signatures/Metanalysis6h24hByGeneSymbolsMerger.Rds")

# when
result <- sut$merge("test/drug_signature/data/signatures/metanalysis/", gene_symbols)
#saveRDS(result, "test/drug_signature/data/signatures/Metanalysis6h24hByGeneSymbolsMerger.Rds")
#write.table(result, "test/drug_signature/data/signatures/Metanalysis6h24hByGeneSymbolsMerger.csv", sep = ";", row.names = F, dec = ",")

# then
test_that("Metanalysis6h24hByGeneSymbolsMergerTest", {
  expect_identical(result, expected)
}
)