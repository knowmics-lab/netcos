library(testthat)

source("modules/config.R")
source("modules/drug_signature/setuper/LINCSExperimentMetaDataSetuper.R")
source("LINCS_export/LINCSExperimentDataRowsLoader.R")
source("LINCS_export/LINCSExportByGeneList.R")

# setup
gene_computed_dir <- "test/LINCS_export/data/computed/"
gene_expected_dir <- "test/LINCS_export/data/expected/"
lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)
lincsExperimentDataRowsLoader <- LINCSExperimentDataRowsLoader$new(drug_signature_cfg$experiments_data_filename)
sut <- LINCSExportByGeneList$new(lincsExperimentDataRowsLoader, gene_computed_dir)

# given
bing_genes <- read.table("data/LINCS-GSE92742/bing_genes.csv", sep = ";", header = T)
experiments_env <- lincsExperimentMetaDataSetuper$setup(c("24"))
block_size <- 3
gene_start_index <- 10
expected22 <- readRDS(paste0(gene_expected_dir, "22.Rds"))
expected23 <- readRDS(paste0(gene_expected_dir, "23.Rds"))
expected25 <- readRDS(paste0(gene_expected_dir, "25.Rds"))

# when
sut$export(bing_genes, experiments_env, gene_start_index, gene_start_index + block_size)

# then

test_that("DrugSignatureDDR1SlowRDSTest", {
  expect_identical(expected22, readRDS(paste0(gene_computed_dir, "22.Rds")))
  expect_identical(expected23, readRDS(paste0(gene_computed_dir, "23.Rds")))
  expect_identical(expected25, readRDS(paste0(gene_computed_dir, "25.Rds")))
}
)