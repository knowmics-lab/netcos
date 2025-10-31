library(testthat)

source("test/test_config.R")
source("modules/config.R")
source("modules/drug_signature/setuper/LINCSExperimentMetaDataSetuper.R")
source("modules/drug_signature/loader/GCTX_LINCSExperimentDataLoader.R")
source("modules/drug_signature/loader/GeneInfoListLoader.R")
source("modules/drug_signature/lmer/LmerLMM.R")
source("modules/drug_signature/lmer/mapper/LmerLMMToDataFrameMapper.R")
source("modules/drug_signature/DrugSignature.R")
source("modules/drug_signature/DrugSignatureByScenarioAsync.R")
source("modules/drug_signature/DrugSignatureByPerturbationTime.R")

# setup
lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)
gctx_LINCSExperimentDataLoader <- GCTX_LINCSExperimentDataLoader$new(drug_signature_cfg$experiments_data_filename)
geneInfoListLoader <- GeneInfoListLoader$new(drug_signature_cfg$gene_info_filename)
lmm <- LmerLMM$new()
lmmToDataFrameMapper <- LmerLMMToDataFrameMapper$new()
drugSignature <- DrugSignature$new(gctx_LINCSExperimentDataLoader, lmm, lmmToDataFrameMapper, test_dir$signatures_base)
drugSignatureByScenarioAsync <- DrugSignatureByScenarioAsync$new(drugSignature)
sut <- DrugSignatureByPerturbationTime$new(lincsExperimentMetaDataSetuper, geneInfoListLoader, drugSignatureByScenarioAsync)

# given
drugs_filter <- c("AM-251", "AM-404", "AM-580", "aminoglutethimide", "aminopurvalanol-a")
gene_symbols <- c("DDR1", "PAX8")
perturbation_time <- "6"
expected1 <- readRDS("test/drug_signature/data/signatures/DDR1_780_6h-expected.Rds")
expected2 <- readRDS("test/drug_signature/data/signatures/PAX8_7849_6h-expected.Rds")

# when
startTime <- Sys.time()
drug_signature_cfg$skip_already_computed_genes <- F
result <- sut$compute(gene_symbols, perturbation_time, drugs_filter)
drug_signature_cfg$skip_already_computed_genes <- T
totalTime <- Sys.time() - startTime
print(sprintf("end overall computation, time: %s %s", totalTime, attr(totalTime, "units")))

# then
test_that("DrugSignatureByPerturbationTimeTest", {
  expect_equal(readRDS("test/drug_signature/data/signatures/DDR1_780_6h.Rds"), expected1)
  expect_equal(readRDS("test/drug_signature/data/signatures/PAX8_7849_6h.Rds"), expected2)
}
)