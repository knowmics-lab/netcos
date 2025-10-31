library(R6)
source("modules/disease_signature/abstract_class/GeneFilter.R")

LogarithmGeneFilter <- R6Class(
  "LogarithmGeneFilter",
  inherit = GeneFilter,
  public = list(
    filter = function(gene_expressions, disease_control_groups) {
      # ci potrebbero essere degli zeri in gene_expressions, + 1
      # serve ad evitare errori nel calcolo del logaritmo
      mean <- rowMeans(log(gene_expressions + 1, 2))
      return(gene_expressions[mean > 6,])
    }
  )
)