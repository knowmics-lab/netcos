library(R6)
library(foreach)
source("modules/drug_signature/abstract_class/DrugSignatureByScenario.R")
source("modules/CoreProcessorInit.R")

DrugSignatureByScenarioAsync <- R6Class(
  "DrugSignatureByScenarioAsync",
  inherit = DrugSignatureByScenario,
  public = list(
    initialize = function(drugSignature) {
      private$drugSignature <- drugSignature
      private$coreProcessorInit <- CoreProcessorInit$new()
    },

    compute = function(geneScenario, experiments_env) {
      totScenarios <- geneScenario$count()
      private$coreProcessorInit$init()
      foreach(i = 1:totScenarios) %dopar% {
        scenario <- geneScenario$get(i)
        private$drugSignature$compute(experiments_env[[scenario$perturbation_time]], scenario$gene_id, scenario$gene_symbol, i)
      }
      return(NA)
    }
  ),
  private = list(
    drugSignature = NA,
    coreProcessorInit = NA
  )
)
