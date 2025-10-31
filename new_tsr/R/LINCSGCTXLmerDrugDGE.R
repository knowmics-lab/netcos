
LINCSGCTXLmerDrugDGE <- R6Class(
  "LINCSGCTXLmerDrugDGE",
  public = list(
    initialize = function(lincs_level3_filename, output_DGE_dir, skip_already_computed_genes = F) {
      lincsMetadataSetuper <- LINCSMetadataSetuper$new()
      lincsGCTXDataLoader <- LINCSGCTXDataLoader$new(lincs_level3_filename)
      drugDGE <- LMMDGEByGene$new(LmerLMM$new(config$LINCSLMMFormula), LINCSLmerLMMToDataFrameMapper$new())
      private$lincsDrugDGE <- LINCSDrugDGE$new(lincsMetadataSetuper, lincsGCTXDataLoader, drugDGE, output_DGE_dir, skip_already_computed_genes)
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
