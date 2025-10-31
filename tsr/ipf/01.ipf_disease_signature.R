source("ipf/config/ipf_cfg.R")
source("modules/disease_signature/filter/LogarithmGeneFilter.R")
source("modules/disease_signature/OverallDiseaseSignature.R")

genefilter <- LogarithmGeneFilter$new()
overallDiseaseSignature <- OverallDiseaseSignature$new(genefilter)

overallDiseaseSignature$compute(
  ipf_cfg$rna_seq_filename,
  ipf_cfg$rna_seq_samples_groups_map,
  ipf_cfg$name,
  ipf_cfg$signature_filename
)