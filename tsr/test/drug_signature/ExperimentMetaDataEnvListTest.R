library(testthat)
source("modules/config.R")
source("modules/drug_signature/setuper/LINCSExperimentMetaDataSetuper.R")
source("modules/drug_signature/ExperimentMetaDataEnvList.R")

# setup
lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)
sut <- ExperimentMetaDataEnvList$new(lincsExperimentMetaDataSetuper)

# given
perturbation_times <- c("6", "24")
drugs_filter <- c("BMS-777607", "BMY-14802")

# when
result <- sut$get(perturbation_times, drugs_filter)

# then
test_that("ExperimentMetaDataEnvListTest", {
  expect_equal(length(result), 2)
  expect_equal(length(result[["6"]]$experiments_meta_data$pert_id), 609)
  expect_equal(length(result[["24"]]$experiments_meta_data$pert_id), 630)
  expect_identical((result[["6"]]$experiments_meta_data$pert_id[7]), "BRD-A15435692")
  expect_identical(result[["6"]]$experiments_meta_data$inst_id[11], "CPC001_HA1E_6H_X2_B3_DUO52HI53LO:B18")
  expect_identical(result[["24"]]$experiments_meta_data$pert_id[1], "DMSO")
}
)