library(testthat)
source("test/test_config.R")
source("modules/drug_signature/julia/factory/JuliaDrugSignatureFactory.R")

# setup
rds_LINCSExperimentDataLoader <- RDS_LINCSExperimentDataLoader$new("test/drug_signature/data/LINCS_splitted_level3/")
sut <- JuliaDrugSignatureFactory$new(rds_LINCSExperimentDataLoader, test_dir$signatures_base)

# given
gene_id <- "2101"
gene_symbol <- "ESRRA"
computation_number <- 1
experiments_env <- lincsExperimentMetaDataSetuper$setup("24")
experiments_env$experiments_meta_data <- experiments_env$experiments_meta_data[3:50,]
expected <- readRDS("test/drug_signature/data/JuliaDrugSignatureExpected.csv.Rds")

# when
drugSignature <- sut$create()
drug_signature_cfg$skip_already_computed_genes <- F
result <- drugSignature$compute(experiments_env, gene_id, gene_symbol, computation_number)
drug_signature_cfg$skip_already_computed_genes <- T

# then
test_that("JuliaDrugSignatureFactoryTest", {
  expect_identical(result, expected)
}
)