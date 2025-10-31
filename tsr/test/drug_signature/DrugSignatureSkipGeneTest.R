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
sut <- DrugSignature$new(rds_LINCSExperimentDataLoader, lmm, lmmToDataFrameMapper, test_dir$genes_signature)

# given
gene_id <- "2597"
gene_symbol <- "GAPDH"
computation_number <- 1
experiments_env <- lincsExperimentMetaDataSetuper$setup("24")

# when
result <- sut$compute(experiments_env, gene_id, gene_symbol, computation_number)
#write.table(result, "test/drug_signature/data/drug_signature_by_gene.csv", sep = ";", row.names = F)
#saveRDS(result, "test/drug_signature/data/drug_signature_by_gene.Rds")

# then
test_that("DrugByGeneSignatureTest", {
  expect_equal(result, NA)
}
)