library(R6)
library(rhdf5)
source("modules/config.R")
source("modules/TSRLogger.R")
source("modules/obj_is_na.R")

DrugSignatureExportExperimentData <- R6Class(
  "DrugSignatureExportExperimentData",
  public = list(
    initialize = function(experimentDataLoader, export_dir = NA) {
      if (!"ExperimentDataLoader" %in% class(experimentDataLoader))
        stop("the experimentDataLoader instance must by of type ExperimentDataLoader")
      private$experimentDataLoader <- experimentDataLoader
      private$tsrLogger <- TSRLogger$new()
      if (obj_is_na(export_dir)) {
        private$export_dir <- paste0(drug_signature_cfg$signature_base_dir, drug_signature_cfg$genes_export_dir)
      }else {
        private$export_dir <- export_dir
      }
    },
    compute = function(experiments_env, gene_id, gene_symbol, computation_number) {
      pert_time_hours <- experiments_env$experiments_meta_data$pert_time[1]
      private$tsrLogger$log(sprintf("start export gene %s, purtubation time hours: %s, number: %s", gene_symbol, pert_time_hours, computation_number))
      startTime <- Sys.time()
      gene_expressions <- private$experimentDataLoader$load(gene_id, experiments_env, computation_number)
      experiments_data <- experiments_env$experiments_meta_data
      experiments_data$gene_expression <- gene_expressions
      filename <- paste0(private$export_dir, gene_symbol, "_", gene_id, "_", pert_time_hours, "h.csv")
      write.table(experiments_data, filename, row.names = F, quote = F)
      totalTime <- Sys.time() - startTime
      private$tsrLogger$log(sprintf("end computation gene %s, purtubation time hours: %s, number: %s, time: %s %s", gene_symbol, pert_time_hours, computation_number, totalTime, attr(totalTime, "units")))
      return(NA)
    }
  ),
  private = list(
    experimentDataLoader = NA,
    tsrLogger = NA,
    export_dir = NA
  )
)
