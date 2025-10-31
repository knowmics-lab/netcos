library(R6)

source("modules/TSRLogger.R")
source("modules/config.R")
source("modules/drug_signature/ExperimentMetaDataEnvList.R")
source("modules/drug_signature/GeneScenario.R")

DrugSignatureByPerturbationTime <- R6Class(
  "DrugSignatureByPerturbationTime",
  public = list(
    initialize = function(experimentMetaDataSetuper, geneInfoListLoader, drugSignatureByScenario) {
      if (!"DrugSignatureByScenario" %in% class(drugSignatureByScenario))
        stop("the drugSignatureByScenario class type must be a subclass of DrugSignatureByScenario")
      private$experimentMetaDataEnvList <- ExperimentMetaDataEnvList$new(experimentMetaDataSetuper)
      private$geneScenario <- GeneScenario$new(geneInfoListLoader)
      private$drugSignatureByScenario <- drugSignatureByScenario
      private$tsrLogger <- TSRLogger$new()
    },

    compute = function(gene_symbols, perturbation_times, drugs_filter = NA) {
      startTime <- Sys.time()
      perturbation_times_label <- paste0(paste(perturbation_times, collapse = "h, "), "h")
      private$tsrLogger$log(sprintf("start drug signatures computation by scenario, pertubation times: %s", perturbation_times_label))
      private$geneScenario$init(gene_symbols, perturbation_times)
      experiments_env <- private$experimentMetaDataEnvList$get(perturbation_times, drugs_filter, drug_signature_cfg$cell_name)
      private$drugSignatureByScenario$compute(private$geneScenario, experiments_env)
      totalTime <- Sys.time() - startTime
      private$tsrLogger$log(sprintf("end drug by gene signatures computation by pertubation times: %s, computation time: %s %s", perturbation_times_label, totalTime, attr(totalTime, "units")))
      return(NA)
    }
  ),
  private = list(
    experimentMetaDataEnvList = NA,
    geneScenario = NA,
    drugSignatureByScenario = NA,
    tsrLogger = NA
  )
)
