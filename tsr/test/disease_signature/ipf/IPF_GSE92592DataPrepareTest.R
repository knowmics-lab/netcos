library(testthat)
source("modules/config.R")
source("modules/disease_signature/ipf/IPF_GSE92592DataPrepare.R")

# setup
sut <- IPF_GSE92592DataPrepare$new()

# given
expected <- readRDS("test/disease_signature/ipf/data/prepared_gse_92592_data_expected.Rds")

# when
result <- sut$prepare(ipf_disease_signature_cfg$gse_92592_filename, ipf_disease_signature_cfg$gse_92592_annot_filename, ipf_disease_signature_cfg$samples_groups_map)
#saveRDS(result, file = "test/disease_signature/ipf/data/prepared_gse_92592_data_expected.Rds")

test_that("IPF_GSE92592DataPrepareTest", {
  expect_identical(result, expected)
}
)