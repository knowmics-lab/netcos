source("als_NYGC/config/als_NYGC_cfg.R")
source("modules/disease_signature/MostSignificantGenes.R")

mostSignificantGenes <- MostSignificantGenes$new()

mostSignificantGenes$compute(
  als_NYGC_cfg$name,
  gsub(".csv", ".Rds", als_NYGC_cfg$signature_filename),
  als_NYGC_cfg$signature_bin_chen_most_significant_genes_filename,
  als_NYGC_cfg$signature_150_most_significant_genes_filename,
  als_NYGC_cfg$signature_150_most_significant_landmark_genes_filename,
  0.7
)

