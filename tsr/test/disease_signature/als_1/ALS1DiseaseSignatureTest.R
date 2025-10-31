library(testthat)

source("modules/disease_signature/als_1/ALS_GSE3307_GPL97DataPrepare.R")
source("modules/disease_signature/DiseaseSignature.R")

# setup
sut <- DiseaseSignature$new()
als_GSE3307_GPL97DataPrepare <- ALS_GSE3307_GPL97DataPrepare$new()

# given
# 1 = Normal, 0 Als
samples_groups_map <- paste0("XXX11000000000XXXXXX0XXXXX1X1XXXXXXXXXXXXXXXX1XXXX",
                             "XXXXXXXXXXXXXXXXXXXXXXX1XXXXXXXXXXXXXXXXXXXXXXXXXX",
                             "XXXXXX1XX11111X111111")
gene_experiments_data <- als_GSE3307_GPL97DataPrepare$prepare("test/disease_signature/als_1/data/ALS-1-GSE3307/GSE3307-GPL97_series_matrix_test.txt", "test/disease_signature/als_1/data/ALS-1-GSE3307/temp-annot", samples_groups_map)
expected <- readRDS(file = "test/disease_signature/als_1/data/als_1_disease_signature_expected.Rds")

# when
result <- sut$compute(gene_experiments_data)
result[3]<-NULL
#saveRDS(result, "test/disease_signature/als_1/data/als_1_disease_signature_expected.Rds")
#write.table(result, file = "test/disease_signature/als_1/data/als_1_disease_signature_expected.csv", sep = ";", row.names = F)

# then
test_that("ALS1DiseaseSignatureTest", {
  for (i in 1:dim(expected)[2]) {
    expect_equal(result[, i], expected[, i])
  }
}
)