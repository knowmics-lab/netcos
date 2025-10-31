library(testthat)
source("modules/drug_signature/Metanalysis6h24hDrugSignature.R")

# setup
sut <- Metanalysis6h24hDrugSignature$new()

# given
drug_by_gene_signatures_6h <- readRDS("test/drug_signature/data/ma_LINCS_drug_by_gene_signatures_6h.Rds")
drug_by_gene_signatures_24h <- readRDS("test/drug_signature/data/ma_LINCS_drug_by_gene_signatures_24h.Rds")
drug_by_gene_signatures_6h_24h <- merge(x = drug_by_gene_signatures_6h, y = drug_by_gene_signatures_24h, by = c("gene", "drug"))
colnames(drug_by_gene_signatures_6h_24h)[3:5] <- c("DE_log2_FC_6h", "std.error_6h", "t.value_6h")
colnames(drug_by_gene_signatures_6h_24h)[6:8] <- c("DE_log2_FC_24h", "std.error_24h", "t.value_24h")
expected <- readRDS("test/drug_signature/data/meta_analysis_drug_by_gene_signature_6h_24h.Rds")

# when
result <- sut$compute(drug_by_gene_signatures_6h_24h[1,])
#saveRDS(result,"test/drug_signature/data/meta_analysis_drug_by_gene_signature_6h_24h.Rds")
#write.table(result, "test/drug_signature/data/meta_analysis_drug_by_gene_signature_6h_24h.csv", sep = ";", row.names = F)

# then
test_that("MetanalysisDrugSignature6h24hTest", {
  expect_identical(result, expected)
}
)