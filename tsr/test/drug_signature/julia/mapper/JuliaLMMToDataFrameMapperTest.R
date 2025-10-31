library(testthat)

source("modules/config.R")
source("modules/drug_signature/setuper/LINCSExperimentMetaDataSetuper.R")
source("modules/drug_signature/loader/RDS_LINCSExperimentDataLoader.R")
source("modules/drug_signature/julia/mapper/JuliaLMMToDataFrameMapper.R")

#setup
lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)
rds_LINCSExperimentDataLoader <- RDS_LINCSExperimentDataLoader$new("test/drug_signature/data/LINCS_splitted_level3/")
sut <- JuliaLMMToDataFrameMapper$new()

# given
gene_id <- "2101"
gene_symbol <- "ESRRA"
computation_number <- 1
experiments_env <- lincsExperimentMetaDataSetuper$setup("24")
experiments_env$experiments_meta_data <- experiments_env$experiments_meta_data[3:50,]
gene_expressions <- rds_LINCSExperimentDataLoader$load(gene_id, experiments_env, 1)
experiments_meta_data <- experiments_env$experiments_meta_data
experiments_meta_data$gene_expressions <- gene_expressions
LMM_output <- julia_call("fit", julia_eval("LinearMixedModel"), gene_expressions ~ pert_iname + (1 | cell_id) + (1 | rna_plate), experiments_meta_data, show_value = F)
expected <- readRDS("test/drug_signature/data/JuliaObjectLMMToDataFrameMapperExpected.Rds")

# when
result <- sut$map(LMM_output, gene_symbol)
#write.table(result, "test/drug_signature/data/JuliaObjectLMMToDataFrameMapperExpected.csv", sep = ";", row.names = F)
#saveRDS(result, "test/drug_signature/data/JuliaObjectLMMToDataFrameMapperExpected.Rds")

# then
test_that("JuliaObjectLMMToDataFrameMapperTest", {
  expect_identical(result, expected)
}
)