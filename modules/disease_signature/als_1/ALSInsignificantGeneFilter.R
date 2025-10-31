library(R6)
ALSInsignificantGeneFilter <- R6Class(
  "ALSInsignificantGeneFilter",
  public = list(
    filter = function(gene_expressions) {
      # ci potrebbero essere degli zeri in gene_expressions, + 1
      # serve ad evitare errori nel calcolo del logaritmo
      gene_expressions$mean <- rowMeans(log(gene_expressions + 1, 2))
      gene_expressions <- subset(gene_expressions, mean > 6)
      gene_expressions$mean <- NULL
      return(gene_expressions)
    }
  )
)