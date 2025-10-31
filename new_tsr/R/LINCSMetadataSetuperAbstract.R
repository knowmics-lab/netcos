
LINCSMetadataSetuperAbstract <- R6Class(
  "LINCSMetadataSetuperAbstract",
  public = list(
    setup = function(perturbation_time, drugs_filter = NA) {
      stop("I'm an abstract method, implement me")
    }
  )
)