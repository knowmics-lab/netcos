library(R6)

GeneSymbolsLoader <- R6Class(
  "GeneSymbolsLoader",
  public = list(
    load = function(gene_symbols_filename) {
      gene_symbols <- read.table(gene_symbols_filename, sep = ";", header = T)
      return(gene_symbols$gene)
    }
  )
)