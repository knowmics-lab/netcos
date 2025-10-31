library(R6)
source("modules/drug_signature/abstract_class/DrugSignatureByScenario.R")

DrugSignatureByScenarioSync <- R6Class(
  "DrugSignatureByScenarioSync",
  inherit = DrugSignatureByScenario,
  public = list(
    initialize = function(drugSignature) {
      private$drugSignature <- drugSignature
    },

    compute = function(geneScenario, experiments_env) {
      totScenarios <- geneScenario$count()
      for (i in 1:totScenarios) {
        scenario <- geneScenario$get(i)
        private$drugSignature$compute(experiments_env[[scenario$perturbation_time]], scenario$gene_id, scenario$gene_symbol, i)
      }
      return(NA)
    }
  ),
  private = list(
    drugSignature = NA
  )
)
