source("modules/config.R")
source("modules/drug_signature/merger/Metanalysis6h24hMerger.R")

# setup
metanalysis6h24hMerger <- Metanalysis6h24hMerger$new()

gene_symbols_filename <- ipf_disease_signature_cfg$ipf_150_most_significant_landmark_genes_filename
metanalysis_filename <-  paste0(drug_signature_cfg$signature_base_dir, ipf_disease_signature_cfg$metanalysis_150_MS_landmark_base_filename, ".Rds")

no_output <- metanalysis6h24hMerger$merge(gene_symbols_filename, drug_signature_cfg$metanalysis_6h_24h_full_path_dir, metanalysis_filename)
