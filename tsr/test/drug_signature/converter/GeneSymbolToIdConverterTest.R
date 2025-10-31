library(testthat)

source("modules/drug_signature/converter/GeneSymbolToIdConverter.R")
source("modules/drug_signature/loader/GeneInfoListLoader.R")

# setup
gene_info_list_filename <- "test/drug_signature/data/GSE92742_Broad_LINCS_gene_info.txt"
geneInfoListLoader <- GeneInfoListLoader$new(gene_info_list_filename)
sut <- GeneSymbolToIdConverter$new(geneInfoListLoader)

# given
gene_symbol_list1 <- c("EPHB3", "ESRRA", "TRADD", "PRPF8", "CAPNS1", "CFL1", "NARS")
expected1 <- c("2049", "2101", "8717", "10594", "826", "1072", "4677")
gene_symbol_list2 <- c("ESRRA", "EPHB3", "TRADD", "PRPF8", "CFL1", "NARS", "CAPNS1")
expected2 <- c("2101", "2049", "8717", "10594", "1072", "4677", "826")

# when
result1 <- sut$convert(gene_symbol_list1)
result2 <- sut$convert(gene_symbol_list2)

# then
test_that("GeneSymbolToIdConverterTest", {
  expect_identical(result1, expected1)
  expect_identical(result2, expected2)
}
)