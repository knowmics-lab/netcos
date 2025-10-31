library(testthat)
source("modules/drug_signature/setuper/LINCSExperimentMetaDataSetuper.R")
source("modules/config.R")

# setup
sut <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)

# given
selected_pert_time <- "6"
expected <- readRDS("test/drug_signature/data/setup_LINCS_experiments_meta_data_6h.Rds")

# when
result <- sut$setup(selected_pert_time)
#saveRDS(result$experiments_meta_data,"test/drug_signature/data/setup_LINCS_experiments_meta_data_6h.Rds")
#write.table(result$experiments_meta_data, "test/drug_signature/data/setup_LINCS_experiments_meta_data_6h.csv", sep = ";",row.names = FALSE)

# then
test_that("LINCSExperimentMetaDataSetuper6hTest", {
  expect_identical(result$experiments_meta_data, expected)
}
)