source("modules/config.R")

als_NYGC_cfg <- list()
als_NYGC_cfg$name <- "als"
als_NYGC_cfg$rna_seq_metadata_filename <- "data/ALS-NYGC-GSE153960/gse153960_metadata.csv"
als_NYGC_cfg$rna_seq_data_filename <- "data/ALS-NYGC-GSE153960/GSE153960_raw_counts_GRCh38.p13_NCBI.tsv.gz"
als_NYGC_cfg$tissue_statuses_to_be_tested <- c("Non-Neurological Control", "ALS Spectrum MND")
als_NYGC_cfg$tissue_statuses_map <- c("Control", "als")
als_NYGC_cfg$formula <- ~tissue_status + (1 | tissue)

als_NYGC_cfg$signature_base_dir <- "output/disease_signature/als_NYGC/"
als_NYGC_cfg$signature_filename <- paste0(als_NYGC_cfg$signature_base_dir, "als_NYGC_signature.csv")

als_NYGC_cfg$signature_most_significant_genes_base_dir <- paste0(als_NYGC_cfg$signature_base_dir, "most_significant_genes/")
als_NYGC_cfg$signature_150_most_significant_genes_filename <- paste0(als_NYGC_cfg$signature_most_significant_genes_base_dir, "als_NYGC_150_most_significant_genes.csv")
als_NYGC_cfg$signature_150_most_significant_landmark_genes_filename <- paste0(als_NYGC_cfg$signature_most_significant_genes_base_dir, "als_NYGC_150_most_significant_landmark_genes.csv")
als_NYGC_cfg$signature_bin_chen_most_significant_genes_filename <- paste0(als_NYGC_cfg$signature_most_significant_genes_base_dir, "als_NYGC_bin_chen_most_significant_genes.csv")


als_NYGC_cfg$connectivity_score_base_dir <- "output/connectivity_score/als_NYGC/"
als_NYGC_cfg$connectivity_score_150_MS_genes_filename <- paste0(als_NYGC_cfg$connectivity_score_base_dir, "als_NYGC_connectivity_score_150_MS_genes.csv")
als_NYGC_cfg$connectivity_score_150_MS_landmark_genes_filename <- paste0(als_NYGC_cfg$connectivity_score_base_dir, "als_NYGC_connectivity_score_150_MS_landmark_genes.csv")
als_NYGC_cfg$connectivity_score_bin_chen_genes_filename <- paste0(als_NYGC_cfg$connectivity_score_base_dir, "ipf_connectivity_score_bin_chen_genes.csv")

als_NYGC_cfg$connectivity_score_PGx_150_MS_genes_filename <- paste0(als_NYGC_cfg$connectivity_score_base_dir, "als_NYGC_connectivity_score_PGx_150_MS_genes.csv")
als_NYGC_cfg$connectivity_score_PGx_bin_chen_genes_filename <- paste0(als_NYGC_cfg$connectivity_score_base_dir, "ipf_connectivity_score_PGx_bin_chen_genes.csv")
