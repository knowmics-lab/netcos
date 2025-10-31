source("modules/config.R")
source("modules/drug_signature/merger/Metanalysis6h24hMerger.R")

# setup
metanalysis6h24hMerger <- Metanalysis6h24hMerger$new()

gene_symbols_filename <- "data/colon_tumor/colon_tumor_150_most_significant_landmark_genes.csv"
metanalysis_filename <- "output/drug_signature/LINCS/colon_tumor_150_most_significant_landmark_genes_metanalysis_6h_24h.Rds"

no_output <- metanalysis6h24hMerger$merge(gene_symbols_filename, drug_signature_cfg$metanalysis_6h_24h_full_path_dir, metanalysis_filename)
