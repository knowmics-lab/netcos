library(R6)
source("modules/config.R")
source("modules/obj_is_na.R")

LINCSExporterByGeneList <- R6Class(
  "LINCSExporterByGeneList",
  public = list(
    initialize = function(lincsExperimentDataRowsLoader, lincs_splitted_level3_dir = NA) {
      private$lincsExperimentDataRowsLoader <- lincsExperimentDataRowsLoader
      if (obj_is_na(lincs_splitted_level3_dir))
        private$lincs_splitted_level3_dir <- lincs_splitted_level3_dir
      else
        private$lincs_splitted_level3_dir <- lincs_splitted_level3_dir
    },
    export = function(genes, experiments_env, start_idx, end_idx) {
      startTime <- Sys.time()
      print(sprintf("start export gene index: %s, group size: %s, time: %s", start_idx, end_idx - start_idx + 1, startTime))
      selected_genes <- as.character(genes[start_idx:end_idx, "gene_id"])
      experiment_data_rows <- private$lincsExperimentDataRowsLoader$load(selected_genes)
      experiment_data_rows <- experiment_data_rows[, experiments_env$experiments_meta_data$inst_id]
      experiment_data_rows_colnames <- colnames(experiment_data_rows)
      for (i in 1:dim(experiment_data_rows)[1]) {
        filename <- paste0(private$lincs_splitted_level3_dir, selected_genes[i], ".Rds")
        experiment_data <- matrix(experiment_data_rows[i,], nrow = 1)
        colnames(experiment_data) <- experiment_data_rows_colnames
        saveRDS(experiment_data, filename)
      }
      totalTime <- Sys.time() - startTime
      print(sprintf("end export gene index: %s, time: %s %s", start_idx, totalTime, attr(totalTime, "units")))
    }
  ),
  private = list(
    lincsExperimentDataRowsLoader = NA,
    lincs_splitted_level3_dir = NA
  )
)