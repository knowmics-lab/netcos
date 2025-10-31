source("als_2/config/als_2_signature_cfg.R")
source("modules/disease_signature/MostSignificantGenes.R")

mostSignificantGenes <- MostSignificantGenes$new()

mostSignificantGenes$compute(
  gsub(".csv", ".Rds", als_2_signature_cfg$signature_filename),
  als_2_signature_cfg$bin_chen_most_significant_genes_filename,
  als_2_signature_cfg$a150_most_significant_genes_filename,
  als_2_signature_cfg$a150_most_significant_landmark_genes_filename
)

