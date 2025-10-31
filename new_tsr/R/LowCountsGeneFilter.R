
LowCountsGeneFilter <- R6Class(
  "LowCountsGeneFilter",
  inherit = GeneFilterAbstract,
  public = list(
    filter = function(gene_expressions, samples) {
      # pre-filter low count genes
      # keep genes with at least N counts > 10, where N = size of smallest group
      keep <- rowSums(gene_expressions >= 10) >= min(table(samples))
      gene_expressions <- gene_expressions[keep,]
      return(gene_expressions)
    }
  )
)