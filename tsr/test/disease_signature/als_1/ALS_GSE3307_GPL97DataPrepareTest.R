library(testthat)
source("modules/disease_signature/als_1/ALS_GSE3307_GPL97DataPrepare.R")

# setup
sut <- ALS_GSE3307_GPL97DataPrepare$new()

# given
samples_groups_map <- paste0("1XX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                             "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                             "XXXXXXXXXXXXXXXXXXXXX")
gse_3307_GPL97_filename <- "test/disease_signature/als_1/data/ALS-1-GSE3307/GSE3307-GPL97_series_matrix_test.txt"
gse_3307_GPL97_temp_annot_dir <- "test/disease_signature/als_1/data/ALS-1-GSE3307/temp-annot"
experiments_count <- 10
sample_count <- 2

# when
result <- sut$prepare(gse_3307_GPL97_filename, gse_3307_GPL97_temp_annot_dir, samples_groups_map)

# then
gene_expressions_matrix <- as.matrix(result$gene_expressions)

test_that("ALS_GSE3307_GPL97DataPrepareTest", {
  expect_equal(length(result), 2)

  expect_equal(dim(gene_expressions_matrix), c(experiments_count, sample_count))
  expect_identical(colnames(gene_expressions_matrix)[1:2], c("GSM120591", "GSM120719"))
  expect_identical(rownames(gene_expressions_matrix)[1:2], c("ARF3", "CAPNS1"))
  expect_identical(gene_expressions_matrix[8][1], 1023)

  expect_equal(length(result$disease_control_groups), 2)
  expect_identical(levels(result$disease_control_groups), c("Control", "ALS"))
}
)