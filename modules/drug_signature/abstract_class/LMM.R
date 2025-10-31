library(R6)

LMM <- R6Class(
  "LMM",
  public = list(
    initialize = function(name) {
      private$name <- name
    },
    compute = function(experiments_env, gene_expressions) {
      stop("I'm an abstract method, implement me")
    },
    algorithmName = function() {
      return(private$name)
    }
  ),
  private = list(
    name = NA
  )
)