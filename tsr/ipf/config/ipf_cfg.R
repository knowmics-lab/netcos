source("modules/config.R")

ipf_cfg <- list()
ipf_cfg$name <- "ipf"
ipf_cfg$rna_seq_filename <- "data/IPF-GSE92592/GSE92592_raw_counts_GRCh38.p13_NCBI.tsv"
ipf_cfg$rna_seq_samples_groups_map <- "000000000000000000001111111111111111111"

ipf_cfg$signature_base_dir <- "output/disease_signature/ipf/"
ipf_cfg$signature_filename <- paste0(ipf_cfg$signature_base_dir, "ipf_signature.csv")

ipf_cfg$signature_most_significant_genes_base_dir <- paste0(ipf_cfg$signature_base_dir, "most_significant_genes/")
ipf_cfg$signature_150_most_significant_genes_filename <- paste0(ipf_cfg$signature_most_significant_genes_base_dir, "ipf_150_most_significant_genes.csv")
ipf_cfg$signature_150_most_significant_landmark_genes_filename <- paste0(ipf_cfg$signature_most_significant_genes_base_dir, "ipf_150_most_significant_landmark_genes.csv")
ipf_cfg$signature_bin_chen_most_significant_genes_filename <- paste0(ipf_cfg$signature_most_significant_genes_base_dir, "ipf_bin_chen_most_significant_genes.csv")

ipf_cfg$connectivity_score_base_dir <- "output/connectivity_score/ipf/"
ipf_cfg$connectivity_score_150_MS_genes_filename <- paste0(ipf_cfg$connectivity_score_base_dir, "ipf_connectivity_score_150_MS_genes.csv")
ipf_cfg$connectivity_score_150_MS_landmark_genes_filename <- paste0(ipf_cfg$connectivity_score_base_dir, "ipf_connectivity_score_150_MS_landmark_genes.csv")
ipf_cfg$connectivity_score_bin_chen_genes_filename <- paste0(ipf_cfg$connectivity_score_base_dir, "ipf_connectivity_score_bin_chen_genes.csv")

