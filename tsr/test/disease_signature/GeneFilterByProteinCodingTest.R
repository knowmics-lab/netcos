library(testthat)
source("modules/disease_signature/GeneFilterByProteinCoding.R")

# setup
sut <- GeneFilterByProteinCoding$new()

# given
gene_expressions <- data.frame(
  esperiment_code = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
)
rownames(gene_expressions) <- c("DDR1", "PAX8", "ABCA1", "ABCA2", "OR4C49P", "ABCB7", "ABCF1", "MIR2117HG", "ABL1", "CYCSP42")
expected <- c("DDR1", "PAX8", "ABCA1", "ABCA2", "ABCB7", "ABCF1", "ABL1")

# when
result <- sut$filter(gene_expressions)

# then
test_that("disease_signature_test", {
  expect_identical(sort(rownames(result)), sort(expected))
}
)