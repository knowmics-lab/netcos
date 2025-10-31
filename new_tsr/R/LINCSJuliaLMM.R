
LINCSJuliaLMM <- R6Class(
  "LINCSJuliaLMM",
  inherit = LMMAbstract,
  public = list(
    compute = function(rna_data) {
      return(
        julia_call("fit", julia_eval("LinearMixedModel"), config$LINCSLMMFormula, rna_data, REML = T)
      )
    }
  )
)
