source("modules/config.R")
source("draft/export/LINCSExperimentMetaDataSetuperMH.R")
source("draft/export/export_group.R")

bing_genes <- read.table("data/LINCS-GSE92742/bing_genes.csv", sep = ";", header = T)
lincsExperimentMetaDataSetuperMH <- LINCSExperimentMetaDataSetuperMH$new()
experiments_env <- lincsExperimentMetaDataSetuperMH$setup(drug_signature_cfg$experiments_meta_data_filename, c("6","24"))
block_size <- 40
gene_index <- 8001
tot_genes <- 9000 #length(bing_genes) = 10174
while (gene_index <= tot_genes) {
  ext_sup <- min(gene_index + block_size - 1, tot_genes)
  export_group(bing_genes, experiments_env, gene_index, ext_sup)
  gene_index <- gene_index + block_size
}