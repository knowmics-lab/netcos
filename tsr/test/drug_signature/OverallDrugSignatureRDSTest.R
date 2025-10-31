library(testthat)

source("test/test_config.R")
source("modules/config.R")
source("modules/drug_signature/loader/GeneInfoListLoader.R")
source("modules/drug_signature/setuper/LINCSExperimentMetaDataSetuper.R")
source("modules/drug_signature/loader/RDS_LINCSExperimentDataLoader.R")
source("modules/drug_signature/lmer/LmerLMM.R")
source("modules/drug_signature/lmer/mapper/LmerLMMToDataFrameMapper.R")
source("modules/drug_signature/DrugSignature.R")
source("modules/drug_signature/DrugSignatureByScenarioAsync.R")
source("modules/drug_signature/DrugSignatureByPerturbationTime.R")
source("modules/drug_signature/OverallDrugSignature.R")

# setup
geneInfoListLoader <- GeneInfoListLoader$new(drug_signature_cfg$gene_info_filename)
lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)
rds_LINCSExperimentDataLoader <- RDS_LINCSExperimentDataLoader$new(drug_signature_cfg$lincs_splitted_level3_dir)
lmm <- LmerLMM$new()
lmmToDataFrameMapper <- LmerLMMToDataFrameMapper$new()
drugSignature <- DrugSignature$new(rds_LINCSExperimentDataLoader, lmm, lmmToDataFrameMapper, test_dir$signatures_base)
drugSignatureByScenarioAsync <- DrugSignatureByScenarioAsync$new(drugSignature)
drugSignatureByPerturbationTime <- DrugSignatureByPerturbationTime$new(lincsExperimentMetaDataSetuper, geneInfoListLoader, drugSignatureByScenarioAsync)
sut <- OverallDrugSignature$new(drugSignatureByPerturbationTime)

# given
drugs_filter <- c("AM-251", "AM-404", "AM-580")
gene_symbols <- c("ADORA3", "AFTPH")
perturbation_times <- c("6", "24")
expected1 <- readRDS("test/drug_signature/data/signatures/ADORA3_140_6h-expected.Rds")
expected2 <- readRDS("test/drug_signature/data/signatures/ADORA3_140_24h-expected.Rds")
expected3 <- readRDS("test/drug_signature/data/signatures/AFTPH_54812_6h-expected.Rds")
expected4 <- readRDS("test/drug_signature/data/signatures/AFTPH_54812_24h-expected.Rds")

# when
drug_signature_cfg$skip_already_computed_genes <- F
sut$compute(gene_symbols, drugs_filter, perturbation_times)
drug_signature_cfg$skip_already_computed_genes <- T

# then
test_that("OverallDrugSignatureTest", {
  expect_identical(readRDS("test/drug_signature/data/signatures/ADORA3_140_6h.Rds"), expected1)
  expect_identical(readRDS("test/drug_signature/data/signatures/ADORA3_140_24h.Rds"), expected2)
  expect_identical(readRDS("test/drug_signature/data/signatures/AFTPH_54812_6h.Rds"), expected3)
  expect_identical(readRDS("test/drug_signature/data/signatures/AFTPH_54812_24h.Rds"), expected4)
}
)
