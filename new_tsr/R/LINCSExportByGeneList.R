
LINCSExportByGeneList <- R6Class(
  "LINCSExportByGeneList",
  public = list(
    initialize = function(lincsDataRowLoader) {
      private$lincsDataRowLoader <- lincsDataRowLoader
    },
    export = function(gene_ids, experiments_meta_data, start_idx, end_idx, output_dir) {
      startTime <- Sys.time()
      dgrpLogger$log(sprintf("start export gene index: %s, group size: %s, time: %s", start_idx, end_idx - start_idx + 1, startTime))
      selected_genes <- as.character(gene_ids[start_idx:end_idx])
      experiment_data_rows <- private$lincsDataRowLoader$load(selected_genes)
      experiment_data_rows <- experiment_data_rows[, experiments_meta_data$inst_id]
      experiment_data_rows_colnames <- colnames(experiment_data_rows)
      for (i in 1:dim(experiment_data_rows)[1]) {
        filename <- paste0(output_dir, selected_genes[i], ".Rds")
        experiment_data <- matrix(experiment_data_rows[i,], nrow = 1)
        colnames(experiment_data) <- experiment_data_rows_colnames
        saveRDS(experiment_data, filename)
      }
      totalTime <- Sys.time() - startTime
      dgrpLogger$log(sprintf("end export gene index: %s, time: %s %s", start_idx, totalTime, attr(totalTime, "units")))
    }
  ),
  private = list(
    lincsDataRowLoader = NA
  )
)