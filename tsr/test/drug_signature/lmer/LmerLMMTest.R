library(testthat)

source("modules/config.R")
source("modules/drug_signature/setuper/LINCSExperimentMetaDataSetuper.R")
source("modules/drug_signature/loader/RDS_LINCSExperimentDataLoader.R")
source("modules/drug_signature/lmer/LmerLMM.R")

#setup
lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)
rds_LINCSExperimentDataLoader <- RDS_LINCSExperimentDataLoader$new("test/drug_signature/data/LINCS_splitted_level3/")
sut <- LmerLMM$new()

# given
gene_id <- "2101"
gene_symbol <- "ESRRA"
computation_number <- 1
experiments_env <- lincsExperimentMetaDataSetuper$setup("24")
experiments_env$experiments_meta_data <- experiments_env$experiments_meta_data[3:50,]
gene_expressions <- rds_LINCSExperimentDataLoader$load(gene_id, experiments_env, 1)
expected <- readRDS("test/drug_signature/data/LmerLMMExpected.Rds")

# when
result <- sut$compute(experiments_env, gene_expressions)
result_table <- coef(summary(result))
result_name <- sut$algorithmName()
#write.table(result_table, "test/drug_signature/data/LmerLMMExpected.csv", sep = ";", row.names = F)
#saveRDS(result_table, "test/drug_signature/data/LmerLMMExpected.Rds")

# then
test_that("LmerLMMTest", {
  expect_identical(result_table, expected)
  expect_identical(result_name, "lmer LMM")
}
)