
JuliaLMM <- R6Class(
  "JuliaLMM",
  inherit = LMMAbstract,
  public = list(
    initialize = function(formula) {
      private$formula <- formula
    },
    compute = function(rna_data) {
      return(
        julia_call("fit", julia_eval("LinearMixedModel"), private$formula, rna_data, REML = T)
      )
    }
  ),
  private = list(
    formula = NA
  )
)
