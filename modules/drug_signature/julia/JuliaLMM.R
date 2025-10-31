library(R6)

source("modules/drug_signature/abstract_class/LMM.R")

JuliaLMM <- R6Class(
  "JuliaLMM",
  inherit = LMM,
  public = list(
    initialize = function() {
      super$initialize("Julia LMM")
    },
    compute = function(experiments_env, gene_expressions) {
      experiments_meta_data <- experiments_env$experiments_meta_data
      experiments_meta_data$gene_expressions <- gene_expressions
      #julia_assign("experiments_meta_data", experiments_meta_data)
      #julia_command("LMM_output = fit!(REML=true,LinearMixedModel(@formula(gene_expressions ~ pert_iname + (1 | cell_id) + (1 | rna_plate)), experiments_meta_data));")
      #LMM_output <- julia_eval("LMM_output")
      LMM_output <- julia_call("fit", julia_eval("LinearMixedModel"), gene_expressions ~ pert_iname + (1 | cell_id) + (1 | rna_plate), experiments_meta_data, REML = T)
      return(LMM_output)
    }
  )
)