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
gene_id <- "2101"
gene_symbol <- "ESRRA"
computation_number <- 1
experiments_env <- lincsExperimentMetaDataSetuper$setup("24")
experiments_env$experiments_meta_data <- experiments_env$experiments_meta_data[3:50,]
expected <- readRDS("test/drug_signature/data/drug_signature_by_gene.Rds")

# when
drug_signature_cfg$skip_already_computed_genes <- F
result <- sut$compute(experiments_env, gene_id, gene_symbol, computation_number)
drug_signature_cfg$skip_already_computed_genes <- T
#write.table(result, "test/drug_signature/data/drug_signature_by_gene.csv", sep = ";", row.names = F)
#saveRDS(result, "test/drug_signature/data/drug_signature_by_gene.Rds")

# then
test_that("DrugByGeneSignatureTest", {
  expect_equal(result, expected)
}
)