library(testthat)

source("modules/config.R")
source("modules/drug_signature/exporter/LINCSExporter.R")

# setup
gene_computed_dir <- "test/LINCS_export/data/computed/"
sut <- LINCSExporter$new(gene_computed_dir)

# given
gene_expected_dir <- "test/LINCS_export/data/expected/"
filename <- "test/drug_signature/data/GCTX_genes_to_export.csv"
expected22 <- readRDS(paste0(gene_expected_dir, "22.Rds"))
expected23 <- readRDS(paste0(gene_expected_dir, "23.Rds"))
expected25 <- readRDS(paste0(gene_expected_dir, "25.Rds"))
expected29 <- readRDS(paste0(gene_expected_dir, "29.Rds"))

# when
sut$export(filename, "24")

# then

test_that("DrugSignatureDDR1SlowRDSTest", {
  expect_identical(expected22, readRDS(paste0(gene_computed_dir, "22.Rds")))
  expect_identical(expected23, readRDS(paste0(gene_computed_dir, "23.Rds")))
  expect_identical(expected25, readRDS(paste0(gene_computed_dir, "25.Rds")))
  expect_identical(expected29, readRDS(paste0(gene_computed_dir, "29.Rds")))
}
)