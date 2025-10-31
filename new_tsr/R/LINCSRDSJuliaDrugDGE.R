
LINCSRDSJuliaDrugDGE <- R6Class(
  "LINCSRDSJuliaDrugDGE",
  public = list(
    initialize = function(lincs_splitted_level3_dir, output_DGE_dir, BLAS_num_threads = NA, skip_already_computed_genes = F) {
      lincsMetadataSetuper <- LINCSMetadataSetuper$new()
      lincsRDSDataLoader <- LINCSRDSDataLoader$new(lincs_splitted_level3_dir)
      juliaLMMDGE <- LINCSJuliaDrugDGEFactory$new()$create(BLAS_num_threads)
      private$lincsDrugDGE <- LINCSDrugDGE$new(lincsMetadataSetuper, lincsRDSDataLoader, juliaLMMDGE, output_DGE_dir, skip_already_computed_genes)
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
