library(testthat)
library(lme4)
source("modules/config.R")
source("modules/drug_signature/LINCSExperimentMetaDataSetuper.R")
source("modules/drug_signature/RDS_LINCSExperimentDataLoader.R")
source("modules/drug_signature/mapper/LmerLMMToDataFrameMapper.R")

#setup
lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)
rds_LINCSExperimentDataLoader <- RDS_LINCSExperimentDataLoader$new("test/drug_signature/data/LINCS_splitted_level3/")
sut <- LmerLMMToDataFrameMapper$new()

# given
gene_id <- "2101"
gene_symbol <- "ESRRA"
computation_number <- 1
experiments_env <- lincsExperimentMetaDataSetuper$setup("24")
experiments_env$experiments_meta_data <- experiments_env$experiments_meta_data[3:50,]
gene_expressions <- rds_LINCSExperimentDataLoader$load(gene_id, experiments_env, 1)
LMM_output <- lmer(gene_expressions ~ pert_iname + (1 | cell_id) + (1 | rna_plate), data = experiments_env$experiments_meta_data, control = lmerControl(calc.derivs = FALSE))
expected <- readRDS("test/drug_signature/data/LmerObjectLMMToDataFrameMapperExpected.Rds")

# when
result <- sut$map(LMM_output, gene_symbol)
#write.table(result, "test/drug_signature/data/LmerObjectLMMToDataFrameMapperExpected.csv", sep = ";", row.names = F)
#saveRDS(result, "test/drug_signature/data/LmerObjectLMMToDataFrameMapperExpected.Rds")

# then
test_that("LmerObjectLMMToDataFrameMapperTest", {
  expect_identical(result, expected)
}
)