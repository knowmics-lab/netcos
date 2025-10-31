library(testthat)

source("test/test_config.R")
source("modules/config.R")
source("modules/drug_signature/setuper/LINCSExperimentMetaDataSetuper.R")
source("modules/drug_signature/loader/RDS_LINCSExperimentDataLoader.R")
source("modules/drug_signature/loader/GeneInfoListLoader.R")
source("modules/drug_signature/lmer/LmerLMM.R")
source("modules/drug_signature/lmer/mapper/LmerLMMToDataFrameMapper.R")
source("modules/drug_signature/DrugSignature.R")
source("modules/drug_signature/DrugSignatureByScenarioAsync.R")
source("modules/drug_signature/DrugSignatureByPerturbationTime.R")

# setup
lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)
rds_LINCSExperimentDataLoader <- RDS_LINCSExperimentDataLoader$new(drug_signature_cfg$lincs_splitted_level3_dir)
geneInfoListLoader <- GeneInfoListLoader$new(drug_signature_cfg$gene_info_filename)
lmm <- LmerLMM$new()
lmmToDataFrameMapper <- LmerLMMToDataFrameMapper$new()
drugSignature <- DrugSignature$new(rds_LINCSExperimentDataLoader, lmm, lmmToDataFrameMapper, test_dir$signatures_base)
drugSignatureByScenarioAsync <- DrugSignatureByScenarioAsync$new(drugSignature)

sut <- DrugSignatureByPerturbationTime$new(lincsExperimentMetaDataSetuper, geneInfoListLoader, drugSignatureByScenarioAsync)

# given
drugs_filter <- c("AM-251", "AM-404", "AM-580", "aminoglutethimide", "aminopurvalanol-a")
gene_symbols <- c("DDR1")
perturbation_time <- "24"
expected <- readRDS("test/drug_signature/data/signatures/DDR1_780_24h-expected.Rds")

# when
startTime <- Sys.time()
drug_signature_cfg$skip_already_computed_genes <- F
sut$compute(gene_symbols, perturbation_time, drugs_filter)
drug_signature_cfg$skip_already_computed_genes <- T
totalTime <- Sys.time() - startTime
print(sprintf("end overall computation, time: %s %s", totalTime, attr(totalTime, "units")))

# then
test_that("DrugSignatureByPerturbationTimeRDS3Test", {
  expect_identical(readRDS("test/drug_signature/data/signatures/DDR1_780_24h.Rds"), expected)
}
)