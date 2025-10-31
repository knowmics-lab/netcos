library(R6)

GeneFilter <- R6Class(
  "GeneFilter",
  public = list(
    filter = function(gene_expressions, disease_control_groups) {
      stop("I'm an abstract method, implement me")
    }
  )
)