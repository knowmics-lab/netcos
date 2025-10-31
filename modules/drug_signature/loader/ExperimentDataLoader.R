library(R6)

ExperimentDataLoader <- R6Class(
  "ExperimentDataLoader",
  public = list(
    load = function(gene_id, experiments_env, computation_number) {
      stop("I'm an abstract method, implement me")
    }
  )
)
