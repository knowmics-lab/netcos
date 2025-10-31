source("modules/config.R")

als_2_signature_cfg <- list()
als_2_signature_cfg$name <- "als"
als_2_signature_cfg$rna_seq_filename <- "data/ALS-2-GSE234297/GSE234297_raw_counts_GRCh38.p13_NCBI.tsv.gz"
als_2_signature_cfg$samples_groups_map <- paste0("00000000000000000000000000000000000000000000000000",
                                                 "00000000000000000000000000000000000011111111111111",
                                                 "11111111111111111111111111111111")
als_2_signature_cfg$output_base_dir <- "output/disease_signature/als_2/"
als_2_signature_cfg$signature_dir <- "als_2/"
als_2_signature_cfg$signature_filename <- paste0(disease_signature_cfg$base_dir, als_2_signature_cfg$signature_dir, "als_2_signature.csv")
als_2_signature_cfg$a150_most_significant_genes_filename <- paste0(als_2_signature_cfg$output_base_dir, "als_2_150_most_significant_genes.csv")
als_2_signature_cfg$a150_most_significant_landmark_genes_filename <- paste0(als_2_signature_cfg$output_base_dir, "als_2_150_most_significant_landmark_genes.csv")
als_2_signature_cfg$bin_chen_most_significant_genes_filename <- paste0(als_2_signature_cfg$output_base_dir, "als_2_bin_chen_most_significant_genes.csv")
als_2_signature_cfg$connectivity_score_150_MS_genes_filename <- paste0(als_2_signature_cfg$output_base_dir, "als_2_connectivity_score_150_MS_genes_filename.csv")
als_2_signature_cfg$connectivity_score_150_MS_landmark_genes_filename <- paste0(als_2_signature_cfg$output_base_dir, "als_2_connectivity_score_150_MS_landmark_genes_filename.csv")

