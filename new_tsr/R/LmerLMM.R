
LmerLMM <- R6Class(
  "LmerLMM",
  inherit = LMMAbstract,
  public = list(
    initialize = function(formula) {
      private$formula <- formula
    },
    compute = function(rna_data) {
      return(
        lmer(private$formula, data = rna_data, control = lmerControl(calc.derivs = FALSE))
      )
    }
  ),
  private = list(
    formula = NA
  )
)