library(R6)
source("modules/drug_signature/loader/GeneSymbolsLoader.R")

GeneSymbolsVectorSplitter <- R6Class(
  "GeneSymbolsVectorSplitter",
  public = list(
    split = function(filename, num_chunks) {
      gene_symbols <- private$geneSymbolsLoader$load(filename)
      size <- length(gene_symbols)
      chunk_size <- ceiling(size / num_chunks)
      chunks <- list()
      start_index <- 1
      if (num_chunks > 1) {
        for (i in 1:(num_chunks - 1)) {
          chunks[[i]] <- gene_symbols[start_index:(chunk_size * i)]
          start_index <- chunk_size * i + 1
        }
      }
      chunks[[num_chunks]] <- gene_symbols[start_index:size]
      return(chunks)
    }
  ),
  private = list(
    geneSymbolsLoader = GeneSymbolsLoader$new()
  )
)