
LINCSDrugDGE <- R6Class(
  "LINCSDrugDGE",
  public = list(
    initialize = function(metadataSetuper, geneRNADataLoader, drugDGE, output_DGE_dir, skip_already_computed_genes = F) {
      private$metadataSetuper <- metadataSetuper
      private$geneRNADataLoader <- geneRNADataLoader
      private$output_DGE_dir <- output_DGE_dir
      private$geneIdToSymbolConverter <- GeneIdSymbolConverter$new()
      private$drugDGE <- drugDGE
      private$skip_already_computed_genes <- skip_already_computed_genes
    },

    compute = function(perturbation_times, gene_list, drugs_filter = NA) {
      tot_genes <- length(gene_list)
      tot_perturbation_times <- length(perturbation_times)
      for (t in 1:tot_perturbation_times) {
        dgrpLogger$log(sprintf("start computation by perturbation time: %sh", perturbation_times[t]))
        metadata <- private$metadataSetuper$setup(perturbation_times[t], drugs_filter)
        for (i in 1:tot_genes) {
          gene_symbol <- private$geneIdToSymbolConverter$idToSymbol(gene_list[i])
          filename <- paste0(private$output_DGE_dir, gene_symbol, "_", gene_list[i], "_", perturbation_times[t], "h.Rds")
          if (private$skip_already_computed_genes & file.exists(filename)) {
            dgrpLogger$log(sprintf("skipping gene %s, purtubation time hours: %s, number: %s", gene_symbol, pert_time_hours, computation_number))
            return(NA)
          }
          rna_data <- private$geneRNADataLoader$load(gene_list[i], metadata)
          dge <- private$drugDGE$compute(rna_data, gene_symbol)
          saveRDS(dge, file = filename)
        }
        dgrpLogger$log(sprintf("end computation by perturbation time: %sh", perturbation_times[t]))
      }
      return(NA)
    }
  ),
  private = list(
    metadataSetuper = NA,
    geneRNADataLoader = NA,
    drugDGE = NA,
    geneIdToSymbolConverter = NA,
    output_DGE_dir = NA,
    skip_already_computed_genes = NA
  )
)
