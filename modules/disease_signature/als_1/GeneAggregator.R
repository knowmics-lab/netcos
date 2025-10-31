library(R6)

GeneAggregator <- R6Class(
  "GeneAggregator",
  public = list(
    aggregate_by_gene = function(gene_expressions) {
      gene_expression_mean <- aggregate(. ~ gene, data = gene_expressions, FUN = max)
      rownames(gene_expression_mean) <- gene_expression_mean$gene
      gene_expression_mean$gene <- NULL
      return(gene_expression_mean)
    }
  )
)