
LINCSRDSLmerDrugDGE <- R6Class(
  "LINCSRDSLmerDrugDGE",
  public = list(
    initialize = function(lincs_splitted_level3_dir, output_dge_dir, skip_already_computed_genes = F) {
      lincsMetadataSetuper <- LINCSMetadataSetuper$new()
      lincsRDSDataLoader <- LINCSRDSDataLoader$new(lincs_splitted_level3_dir)
      drugDGE <- LMMDGEByGene$new(LmerLMM$new(config$LINCSLMMFormula), LINCSLmerLMMToDataFrameMapper$new())
      private$lincsDrugDGE <- LINCSDrugDGE$new(lincsMetadataSetuper, lincsRDSDataLoader, drugDGE, output_dge_dir, skip_already_computed_genes)
    },

    compute = function(perturbation_times, gene_list, drugs_filter = NA) {
      private$lincsDrugDGE$compute(perturbation_times, gene_list, drugs_filter)
      return(NA)
    }
  ),
  private = list(
    lincsDrugDGE = NA
  )
)
