source("als_NYGC/config/als_NYGC_cfg.R")
source("modules/disease_signature/filter/LogarithmGeneFilter.R")
source("modules/disease_signature/LMM/OverallDiseaseSignatureLMM.R")

geneFilter <- LogarithmGeneFilter$new()
overallDiseaseSignatureLMM <- OverallDiseaseSignatureLMM$new(geneFilter)

overallDiseaseSignatureLMM$compute(
  als_NYGC_cfg$rna_seq_metadata_filename,
  als_NYGC_cfg$rna_seq_data_filename,
  als_NYGC_cfg$tissue_statuses_to_be_tested,
  als_NYGC_cfg$tissue_statuses_map,
  als_NYGC_cfg$formula,
  als_NYGC_cfg$signature_filename
)