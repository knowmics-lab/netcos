library(R6)
source("modules/config.R")
source("modules/connectivity_score/GeneSelectionByDiseaseGeneCount.R")
GeneSelectionStrategy <- R6Class(
  "GeneSelectionStrategy",
  public = list(
    compute = function(disease_gene_count, gene_selection_strategy) {
      if (gene_selection_strategy == disease_signature_cfg$gene_select_strategy_bin_chen) {
        gene_selection <- data.frame(
          name = "bin_chen",
          count = disease_gene_count
        )
      }else {
        gene_selection <- private$geneSelectionByDiseaseGeneCount$compute(disease_gene_count)
      }
      return(gene_selection)
    }
  ),
  private = list(
    geneSelectionByDiseaseGeneCount = GeneSelectionByDiseaseGeneCount$new()
  )
)