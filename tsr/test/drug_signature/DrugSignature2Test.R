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
sut <- DrugSignature$new(rds_LINCSExperimentDataLoader, lmm, lmmToDataFrameMapper, test_dir$signatures_base)

# given
gene_id <- "780"
gene_symbol <- "DDR1"
computation_number <- 1
drugs_filter <- c("AM-251", "AM-404", "AM-580", "aminoglutethimide", "aminopurvalanol-a")
experiments_env <- lincsExperimentMetaDataSetuper$setup("6", drugs_filter)
expected <- readRDS("test/drug_signature/data/drug_signature_by_gene_2.Rds")

# when
drug_signature_cfg$skip_already_computed_genes <- F
result <- sut$compute(experiments_env, gene_id, gene_symbol, computation_number)
drug_signature_cfg$skip_already_computed_genes <- T
#write.table(result, "test/drug_signature/data/drug_signature_by_gene_2.csv", sep = ";", row.names = F)
#saveRDS(result, "test/drug_signature/data/drug_signature_by_gene_2.Rds")

# then
test_that("DrugByGeneSignatureTest", {
  expect_identical(result, expected)
}
)