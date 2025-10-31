library(testthat)
source("modules/config.R")
source("modules/disease_signature/Illumina_HiSeq_mapper/IlluminaHiSeqMapper.R")
source("modules/disease_signature/DiseaseSignature.R")

# setup
illuminaHiSeqMapper <- IlluminaHiSeqMapper$new()
sut <- DiseaseSignature$new()

# given
gene_experiments_data <- illuminaHiSeqMapper$map(ipf_disease_signature_cfg$gse_92592_filename, ipf_disease_signature_cfg$samples_groups_map,"ipf")
expected <- readRDS(file = "test/disease_signature/ipf/data/ipf_disease_signature_expected.Rds")

# when
result <- sut$compute(gene_experiments_data, F)
#write.table(result, file = "test/disease_signature/ipf/data/ipf_disease_signature_expected.csv", sep = ";", row.names = F)
#saveRDS(result, file = "test/disease_signature/ipf/data/ipf_disease_signature_expected.Rds")

# then
test_that("IPFDiseaseSignatureTest", {
  for (i in 1:dim(expected)[2]) {
    expect_equal(result[, i], expected[, i])
  }
}
)