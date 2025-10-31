
LogarithmGeneFilter <- R6Class(
  "LogarithmGeneFilter",
  inherit = GeneFilterAbstract,
  public = list(
    filter = function(gene_expressions, samples) {
      # ci potrebbero essere degli zeri in gene_expressions,
      #  + 1 serve ad evitare errori nel calcolo del logaritmo
      mean <- rowMeans(log(gene_expressions + 1, 2))
      return(gene_expressions[mean > 6,])
    }
  )
)