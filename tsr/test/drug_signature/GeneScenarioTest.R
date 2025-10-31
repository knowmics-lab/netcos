library(testthat)
source("modules/config.R")
source("modules/drug_signature/GeneScenario.R")
source("modules/drug_signature/loader/GeneInfoListLoader.R")

# setup
gene_ids <- c("2", "  9", "10", "12", "13", "14", "15")
gene_symbols <- c("A2M", "NAT1", "NAT2", "SERPINA3", "AADAC", "AAMP", "AANAT")
perturbation_times <- c("6", "24")
geneInfoListLoader <- GeneInfoListLoader$new(drug_signature_cfg$gene_info_filename)

sut1 <- GeneScenario$new(geneInfoListLoader)
sut1$init(gene_symbols = gene_symbols, perturbation_times = perturbation_times, gene_ids = gene_ids)
sut2 <- GeneScenario$new(geneInfoListLoader)
sut2$init(gene_symbols = gene_symbols, perturbation_times = perturbation_times)

# given
expected1 <- list(gene_id = "2", gene_symbol = "A2M", perturbation_time = "6")
expected2 <- list(gene_id = "15", gene_symbol = "AANAT", perturbation_time = "6")
expected3 <- list(gene_id = "13", gene_symbol = "AADAC", perturbation_time = "6")
expected4 <- list(gene_id = "2", gene_symbol = "A2M", perturbation_time = "24")
expected5 <- list(gene_id = "10", gene_symbol = "NAT2", perturbation_time = "24")
expected6 <- list(gene_id = "9", gene_symbol = "NAT1", perturbation_time = "6")

# when
result1 <- sut1$get(1)
result2 <- sut1$get(7)
result3 <- sut1$get(5)
result4 <- sut1$get(8)
result5 <- sut1$get(10)

result6 <- sut2$get(2)

# then
test_that("GeneScenarioMakerTest", {
  expect_identical(result1, expected1)
  expect_identical(result2, expected2)
  expect_identical(result3, expected3)
  expect_identical(result4, expected4)
  expect_identical(result5, expected5)
  expect_identical(result6, expected6)
}
)