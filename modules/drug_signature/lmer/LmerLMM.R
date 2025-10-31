library(R6)
library(lme4)

source("modules/drug_signature/abstract_class/LMM.R")

LmerLMM <- R6Class(
  "LmerLMM",
  inherit = LMM,
  public = list(
    initialize = function() {
      super$initialize("lmer LMM")
    },
    compute = function(experiments_env, gene_expressions) {
      # (1 | cell_id) +
      LMM_output <- lmer(gene_expressions ~ pert_iname +  (1 | rna_plate), data = experiments_env$experiments_meta_data, control = lmerControl(calc.derivs = FALSE))
      return(LMM_output)
    }
  )
)