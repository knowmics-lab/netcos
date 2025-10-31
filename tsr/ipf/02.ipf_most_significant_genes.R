source("ipf/config/ipf_cfg.R")
source("modules/disease_signature/MostSignificantGenes.R")

mostSignificantGenes <- MostSignificantGenes$new()

mostSignificantGenes$compute(
  ipf_cfg$name,
  gsub(".csv", ".Rds", ipf_cfg$signature_filename),
  ipf_cfg$signature_bin_chen_most_significant_genes_filename,
  ipf_cfg$signature_150_most_significant_genes_filename,
  ipf_cfg$signature_150_most_significant_landmark_genes_filename
)

