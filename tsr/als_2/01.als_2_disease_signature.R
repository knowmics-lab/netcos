source("als_2/config/als_2_signature_cfg.R")
source("modules/disease_signature/OverallDiseaseSignature.R")

# setup
# download genetic db from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE3nnn/GSE3307/matrix/
overallDiseaseSignature <- OverallDiseaseSignature$new()

overallDiseaseSignature$compute(
  als_2_signature_cfg$rna_seq_filename,
  als_2_signature_cfg$samples_groups_map,
  als_2_signature_cfg$name,
  als_2_signature_cfg$signature_filename
)