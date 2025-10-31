library(testthat)

source("test/test_config.R")
source("modules/config.R")
source("modules/drug_signature/lmer/LmerLMM.R")
source("modules/drug_signature/lmer/mapper/LmerLMMToDataFrameMapper.R")
source("modules/drug_signature/DrugSignature.R")
source("modules/drug_signature/setuper/LINCSExperimentMetaDataSetuper.R")
source("modules/drug_signature/loader/RDS_LINCSExperimentDataLoader.R")

# setup
lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)
rds_LINCSExperimentDataLoader <- RDS_LINCSExperimentDataLoader$new("test/drug_signature/data/LINCS_splitted_level3/")
lmm <- LmerLMM$new()
lmmToDataFrameMapper <- LmerLMMToDataFrameMapper$new()
sut <- DrugSignature$new(rds_LINCSExperimentDataLoader, lmm, lmmToDataFrameMapper, test_dir$signatures_base)


# given
drug_ranef_and_sd <- readRDS("test/drug_signature/data/drug_ranef_and_sd.Rds")
PRISM_drugs <- unique(drug_ranef_and_sd$drug_name)
gene_id <- "780"
gene_symbol <- "DDR1"
computation_number <- 1

experiments_env <- lincsExperimentMetaDataSetuper$setup("6", PRISM_drugs)
expected <- readRDS("test/drug_signature/data/drug_by_gene_signature_DDR1_slow.Rds")

# when
drug_signature_cfg$skip_already_computed_genes <- F
result <- sut$compute(experiments_env, gene_id, gene_symbol, computation_number)
drug_signature_cfg$skip_already_computed_genes <- T
#saveRDS(result,"test/drug_signature/data/drug_by_gene_signature_DDR1_slow.Rds")
#write.table(result, "test/drug_signature/data/drug_by_gene_signature_DDR1_slow.csv", sep = ";", row.names = F)

# then
test_that("DrugSignatureDDR1SlowRDSTest", {
  expect_identical(result, expected)
}
)