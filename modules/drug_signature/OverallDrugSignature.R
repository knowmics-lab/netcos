library(R6)

source("modules/config.R")
source("modules/TSRLogger.R")

OverallDrugSignature <- R6Class(
  "OverallDrugSignature",
  public = list(
    initialize = function(drugSignatureByPerturbationTime) {
      private$tsrLogger <- TSRLogger$new()
      private$drugSignatureByPerturbationTime <- drugSignatureByPerturbationTime
    },
    compute = function(gene_symbols, drugs_filter, perturbation_times) {
      startTime <- Sys.time()
      private$tsrLogger$log("start drug signature computation")
      private$drugSignatureByPerturbationTime$compute(gene_symbols, perturbation_times, drugs_filter)
      totalTime <- Sys.time() - startTime
      private$tsrLogger$log(sprintf("end drug signature computation, time: %s %s", totalTime, attr(totalTime, "units")))
      return(NA)
    }
  ),
  private = list(
    drugSignatureByPerturbationTime = NA,
    tsrLogger = NA
  )
)
