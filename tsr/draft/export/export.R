source("modules/config.R")
source("draft/export/LINCSExperimentMetaDataSetuperMH.R")
source("draft/export/load_LINCS_experiment_data_rows.R")

bing_genes <- read.table("data/LINCS-GSE92742/bing_genes.csv", sep = ";", header = T)
selected_genes <- as.character(bing_genes[91:110, "pr_gene_id"])
lincsExperimentMetaDataSetuperMH <- LINCSExperimentMetaDataSetuperMH$new()
experiments_env <- lincsExperimentMetaDataSetuperMH$setup(drug_signature_cfg$experiments_meta_data_filename, c("6", "24"))
experiment_data_rows <- load_LINCS_experiment_data_rows(experiments_env, drug_signature_cfg$experiments_data_filename, selected_genes)
experiment_data_rows_colnames <- colnames(experiment_data_rows)
for (i in 1:dim(experiment_data_rows)[1]) {
  filename <- paste0("data/LINCS-GSE92742/splitted_level3/", selected_genes[i], ".Rds")
  experiment_data <- matrix(experiment_data_rows[i,], nrow = 1)
  colnames(experiment_data) <- experiment_data_rows_colnames
  saveRDS(experiment_data, filename)
}