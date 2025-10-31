library(R6)
source("modules/disease_signature/abstract_class/GeneFilter.R")

LowCountsGeneFilter <- R6Class(
  "LowCountsGeneFilter",
  inherit = GeneFilter,
  public = list(
    filter = function(gene_expressions, disease_control_groups) {
      # pre-filter low count genes
      # keep genes with at least N counts > 10, where N = size of smallest group
      keep <- rowSums(gene_expressions >= 10) >= min(table(disease_control_groups))
      gene_expressions <- gene_expressions[keep,]
      return(gene_expressions)
    }
  )
)