library(R6)

DrugSignatureByScenario <- R6Class(
  "DrugSignatureByScenario",
  public = list(
    compute = function(geneScenario, experiments_env) {
      stop("I'm an abstract method, implement me")
    }
  )
)
