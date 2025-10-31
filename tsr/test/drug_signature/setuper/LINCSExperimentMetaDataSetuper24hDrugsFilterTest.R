library(testthat)
source("modules/drug_signature/setuper/LINCSExperimentMetaDataSetuper.R")
source("modules/config.R")

# setup
sut <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)

# given
selected_pert_time <- "24"
expected <- readRDS("test/drug_signature/data/setup_LINCS_experiments_meta_data_24h_drugs_filter.Rds")
drugs <- c("1-phenylbiguanide", "10-hydroxycamptothecin", "7-nitroindazole", "abiraterone",
           "ABT-737", "ABT-751", "afatinib", "AG-14361", "albendazole", "alfacalcidol",
           "altretamine", "alvespimycin", "alvocidib", "AM-251", "AM-404", "AM-580",
           "aminoglutethimide", "aminopurvalanol-a")

# when
result <- sut$setup(selected_pert_time, drugs)
#saveRDS(result$experiments_meta_data,"test/drug_signature/data/setup_LINCS_experiments_meta_data_24h_drugs_filter.Rds")
#write.table(result$experiments_meta_data, "test/drug_signature/data/setup_LINCS_experiments_meta_data_24h_drugs_filter.csv", sep = ";",row.names = FALSE)

# then
test_that("LINCSExperimentMetaDataSetuper24hDrugsFilterTest", {
  expect_identical(result$experiments_meta_data, expected)
}
)

