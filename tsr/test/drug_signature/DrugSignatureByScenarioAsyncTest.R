library(testthat)

source("test/test_config.R")
source("modules/config.R")
source("modules/drug_signature/setuper/LINCSExperimentMetaDataSetuper.R")
source("modules/drug_signature/ExperimentMetaDataEnvList.R")
source("modules/drug_signature/loader/RDS_LINCSExperimentDataLoader.R")
source("modules/drug_signature/loader/GeneInfoListLoader.R")
source("modules/drug_signature/GeneScenario.R")
source("modules/drug_signature/lmer/LmerLMM.R")
source("modules/drug_signature/lmer/mapper/LmerLMMToDataFrameMapper.R")
source("modules/drug_signature/DrugSignature.R")
source("modules/drug_signature/DrugSignatureByScenarioAsync.R")

# setup
lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)
experimentMetaDataEnvList <- ExperimentMetaDataEnvList$new(lincsExperimentMetaDataSetuper)
rds_LINCSExperimentDataLoader <- RDS_LINCSExperimentDataLoader$new(drug_signature_cfg$lincs_splitted_level3_dir)
geneInfoListLoader <- GeneInfoListLoader$new(drug_signature_cfg$gene_info_filename)
geneScenario <- GeneScenario$new(geneInfoListLoader)
lmm <- LmerLMM$new()
lmmToDataFrameMapper <- LmerLMMToDataFrameMapper$new()
drugSignature <- DrugSignature$new(rds_LINCSExperimentDataLoader, lmm, lmmToDataFrameMapper, test_dir$signatures_base)
sut <- DrugSignatureByScenarioAsync$new(drugSignature)

# given
drugs_filter <- c("AM-251", "AM-404", "AM-580")
gene_symbols <- c("ADORA3", "AFTPH")
perturbation_times <- c("6", "24")
geneScenario$init(gene_symbols, perturbation_times)
experiments_env <- experimentMetaDataEnvList$get(perturbation_times, drugs_filter)
expected1 <- readRDS("test/drug_signature/data/signatures/ADORA3_140_6h-expected.Rds")
expected2 <- readRDS("test/drug_signature/data/signatures/ADORA3_140_24h-expected.Rds")
expected3 <- readRDS("test/drug_signature/data/signatures/AFTPH_54812_6h-expected.Rds")
expected4 <- readRDS("test/drug_signature/data/signatures/AFTPH_54812_24h-expected.Rds")

# when
startTime <- Sys.time()
drug_signature_cfg$skip_already_computed_genes <- F
sut$compute(geneScenario, experiments_env)
drug_signature_cfg$skip_already_computed_genes <- T
totalTime <- Sys.time() - startTime
print(sprintf("end overall computation, time: %s %s", totalTime, attr(totalTime, "units")))

# then
test_that("DrugSignatureByScenarioAsyncTest", {
  expect_equal(readRDS("test/drug_signature/data/signatures/ADORA3_140_6h.Rds"), expected1)
  expect_equal(readRDS("test/drug_signature/data/signatures/ADORA3_140_24h.Rds"), expected2)
  expect_equal(readRDS("test/drug_signature/data/signatures/AFTPH_54812_6h.Rds"), expected3)
  expect_equal(readRDS("test/drug_signature/data/signatures/AFTPH_54812_24h.Rds"), expected4)
}
)