
LMMAbstract <- R6Class(
  "LMMAbstract",
  public = list(
    compute = function(rna_data) {
      stop("I'm an abstract method, implement me")
    }
  )
)