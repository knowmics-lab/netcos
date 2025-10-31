library(testthat)

source("modules/config.R")
source("modules/drug_signature/loader/RDS_LINCSExperimentDataLoader.R")
source("modules/drug_signature/setuper/LINCSExperimentMetaDataSetuper.R")

# setup
sut <- RDS_LINCSExperimentDataLoader$new("test/drug_signature/data/LINCS_splitted_level3/")
lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)

# given
gene_id <- "2101"
computation_number <- 1
experiments_env <- lincsExperimentMetaDataSetuper$setup("24")
experiments_env$experiments_meta_data <- experiments_env$experiments_meta_data[3:50,]
expected <- readRDS("test/drug_signature/data/LINCSExperimentDataLoaderTestExpected.Rds")

# when
result <- sut$load("2101", experiments_env, 24)
#write.table(result, "test/drug_signature/data/LINCSExperimentDataLoaderTestExpected.csv", sep = ";", row.names = F)
#saveRDS(result, "test/drug_signature/data/LINCSExperimentDataLoaderTestExpected.Rds")

# then
test_that("LINCSExperimentDataLoaderTest", {
  expect_identical(result, expected)
}
)