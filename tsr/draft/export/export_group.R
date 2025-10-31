source("modules/config.R")
source("draft/export/LINCSExperimentMetaDataSetuperMH.R")
source("draft/export/load_LINCS_experiment_data_rows.R")

export_group <- function(bing_genes, experiments_env, start_idx, end_idx) {
  startTime <- Sys.time()
  print(sprintf("start export gene index: %s, group size: %s, time: %s", start_idx, end_idx - start_idx + 1, startTime))
  selected_genes <- as.character(bing_genes[start_idx:end_idx, "pr_gene_id"])
  experiment_data_rows <- load_LINCS_experiment_data_rows(experiments_env, drug_signature_cfg$experiments_data_filename, selected_genes)
  experiment_data_rows_colnames <- colnames(experiment_data_rows)
  for (i in 1:dim(experiment_data_rows)[1]) {
    filename <- paste0(drug_signature_cfg$lincs_splitted_level3_dir, selected_genes[i], ".Rds")
    experiment_data <- matrix(experiment_data_rows[i,], nrow = 1)
    colnames(experiment_data) <- experiment_data_rows_colnames
    saveRDS(experiment_data, filename)
  }
  totalTime <- Sys.time() - startTime
  print(sprintf("end export gene index: %s, time: %s %s", start_idx, totalTime, attr(totalTime, "units")))
}