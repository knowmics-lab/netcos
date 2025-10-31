gene_meta_info_cfg <- list()
gene_meta_info_cfg$filename <- "data/gene_meta.Rds"

disease_signature_cfg <- list()
disease_signature_cfg$gene_select_strategy_150_most_significant <- "gene_select_strategy_150_most_significant"
disease_signature_cfg$gene_select_strategy_bin_chen <- "gene_select_strategy_bin_chen"
disease_signature_cfg$Illumina_HiSeq_annot_filename <- "data/Illumina_HiSeq/Human.GRCh38.p13.annot.tsv.gz"
disease_signature_cfg$base_dir <- "tsr/output/disease_signature/"

als_1_disease_signature_cfg <- list()
als_1_disease_signature_cfg$name <- "ALS"
als_1_disease_signature_cfg$gse_3307_GPL97_filename <- "data/ALS-1-GSE3307/GSE3307-GPL97_series_matrix.txt"
als_1_disease_signature_cfg$gse_3307_GPL97_annot_dir <- "data/ALS-1-GSE3307/temp-annot"
als_1_disease_signature_cfg$signature_base_dir <- "tsr/output/disease_signature/als_1/"
als_1_disease_signature_cfg$signature_filename <- "als_signature"
als_1_disease_signature_cfg$als_1_150_most_significant_genes_filename <- paste0(als_1_disease_signature_cfg$signature_base_dir, "als_1_150_most_significant_genes.csv")
als_1_disease_signature_cfg$als_1_150_most_significant_landmark_genes_filename <- paste0(als_1_disease_signature_cfg$signature_base_dir, "als_1_150_most_significant_landmark_genes.csv")
als_1_disease_signature_cfg$als_1_bin_chen_most_significant_genes_filename <- paste0(als_1_disease_signature_cfg$signature_base_dir, "als_bin_chen_most_significant_genes.csv")
# 1 = Control, 0 Als
als_1_disease_signature_cfg$samples_groups_map <- paste0("XXX11000000000XXXXXXXXXXXX1X1XXXXXXXXXXXXXXXX1XXXX",
                                                         "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                                                         "XXXXXX111111111111111")

ipf_disease_signature_cfg <- list()
# 1 = Control, 0 ipf
ipf_disease_signature_cfg$samples_groups_map <- "000000000000000000001111111111111111111"
ipf_disease_signature_cfg$name <- "IPF"
ipf_disease_signature_cfg$gse_92592_filename <- "data/IPF-GSE92592/GSE92592_raw_counts_GRCh38.p13_NCBI.tsv"
ipf_disease_signature_cfg$gse_92592_annot_filename <- "data/IPF-GSE92592/Human.GRCh38.p13.annot.tsv"
ipf_disease_signature_cfg$signature_base_dir <- "tsr/output/disease_signature/ipf/"
ipf_disease_signature_cfg$signature_filename <- "ipf_signature"
ipf_disease_signature_cfg$ipf_150_most_significant_genes_filename <- paste0(ipf_disease_signature_cfg$signature_base_dir, "ipf_150_most_significant_genes.csv")
ipf_disease_signature_cfg$ipf_150_most_significant_landmark_genes_filename <- paste0(ipf_disease_signature_cfg$signature_base_dir, "ipf_150_most_significant_landmark_genes.csv")
ipf_disease_signature_cfg$ipf_bin_chen_most_significant_genes_filename <- paste0(ipf_disease_signature_cfg$signature_base_dir, "ipf_bin_chen_most_significant_genes.csv")
ipf_disease_signature_cfg$drug_signature_genes_by_150_MS_base_filename <- "ipf_drug_signature_150_MS_"
ipf_disease_signature_cfg$drug_signature_genes_by_bin_chen_base_filename <- "ipf_drug_signature_bin_chen_"
ipf_disease_signature_cfg$metanalysis_150_MS_base_filename <- "ipf_metanalysis_150_MS_6h_24h"
ipf_disease_signature_cfg$metanalysis_150_MS_landmark_base_filename <- "ipf_landmark_metanalysis_150_MS_6h_24h"
ipf_disease_signature_cfg$metanalysis_bin_chen_base_filename <- "ipf_metanalysis_bin_chen_6h_24h"
ipf_disease_signature_cfg$connectivity_score_150_MS_filename <- "ipf_connectivity_score_150_MS"
ipf_disease_signature_cfg$connectivity_score_150_MS_landmark_filename <- "ipf_connectivity_score_150_MS_landmark"
ipf_disease_signature_cfg$connectivity_score_bin_chen_filename <- "ipf_connectivity_score_bin_chen"

drug_signature_cfg <- list()
# 65 GB file, download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742 seperately and specify path:
drug_signature_cfg$cell_name <- 'MCF7' # probably wont work with full dataset now...
drug_signature_cfg$lincs_base_dir <- "data/LINCS-GSE92742/"
drug_signature_cfg$drug_list <- "data/drugs.csv"
drug_signature_cfg$gene_info_filename <- paste0(drug_signature_cfg$lincs_base_dir, "GSE92742_Broad_LINCS_gene_info.txt")
drug_signature_cfg$experiments_meta_data_filename <- paste0(drug_signature_cfg$lincs_base_dir, "GSE92742_Broad_LINCS_inst_info.txt")
drug_signature_cfg$experiments_data_filename <- paste0(drug_signature_cfg$lincs_base_dir, "GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx")
drug_signature_cfg$lincs_splitted_level3_dir <- paste0(drug_signature_cfg$lincs_base_dir, "splitted_level3/")
drug_signature_cfg$signature_base_dir <- paste0("tsr/output/drug_signature/LINCS_", drug_signature_cfg$cell_name)
drug_signature_cfg$genes_signature_dir <- "genes/"
drug_signature_cfg$genes_signature_full_path_dir <- paste0(drug_signature_cfg$signature_base_dir, drug_signature_cfg$genes_signature_dir)
drug_signature_cfg$metanalysis_6h_24_dir <- "metanalysis/"
drug_signature_cfg$metanalysis_6h_24h_full_path_dir <- paste0(drug_signature_cfg$signature_base_dir, drug_signature_cfg$metanalysis_6h_24_dir)
drug_signature_cfg$perturbation_times <- c("6", "24")
drug_signature_cfg$gene_symbols_filename <- "data/LINCS-GSE92742/bing_gene_symbols.csv"
drug_signature_cfg$skip_already_computed_genes <- T
drug_signature_cfg$export_block_size <- 100
drug_signature_cfg$bing_genes_filename <- "data/LINCS-GSE92742/bing_genes.csv"

connectivity_score_cfg <- list()
connectivity_score_cfg$cs_base_dir <- "tsr/output/connectivity_score/"
connectivity_score_cfg$cs_filename <- "connectivity_score"

tsr_logger_cfg <- list()
tsr_logger_cfg$log_file_name <- "tsr/log_data/tsr_log_file.txt"

parallel_computation <- list()
parallel_computation$max_cores <- 8
parallel_computation$LMM_disease_signature_max_cores <- 8
parallel_computation$delay_start_clusters <- 40

tsr_julia_cfg <- list()
tsr_julia_cfg$BLAS_num_threads <- 16

if (file.exists("modules/local_config.R")) {
  source("modules/local_config.R")
}
