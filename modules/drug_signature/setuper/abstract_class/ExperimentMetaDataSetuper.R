library(R6)

ExperimentMetaDataSetuper <- R6Class(
  "ExperimentMetaDataSetuper",
  public = list(
    setup = function(perturbation_time, drugs_filter = NA, cell_name = NA) {
      stop("I'm an abstract method, implement me")
    }
  )
)