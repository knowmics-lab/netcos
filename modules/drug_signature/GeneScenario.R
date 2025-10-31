library(R6)
source("modules/obj_is_na.R")
source("modules/drug_signature/converter/GeneSymbolToIdConverter.R")

GeneScenario <- R6Class(
  "GeneScenario",
  public = list(
    initialize = function(geneInfoListLoader) {
      private$geneInfoListLoader <- geneInfoListLoader
    },
    init = function(gene_symbols, perturbation_times, gene_ids = NA) {
      if (obj_is_na(gene_ids)) {
        geneSymbolToIdConverter <- GeneSymbolToIdConverter$new(private$geneInfoListLoader)
        private$gene_ids <- geneSymbolToIdConverter$convert(gene_symbols)
      }else {
        private$gene_ids <- gene_ids
      }
      private$gene_ids_size <- length(private$gene_ids)
      private$gene_symbols <- gene_symbols
      private$perturbation_times <- perturbation_times
    },
    get = function(index) {
      perturbation_time_index <- (index - 1) %/% private$gene_ids_size + 1
      index <- index - private$gene_ids_size * (perturbation_time_index - 1)
      return(list(
        gene_id = private$gene_ids[index],
        gene_symbol = private$gene_symbols[index],
        perturbation_time = private$perturbation_times[perturbation_time_index]
      ))
    },
    count = function() {
      return(private$gene_ids_size * length(private$perturbation_times))
    }
  ),
  private = list(
    gene_ids = NA,
    gene_symbols = NA,
    perturbation_times = NA,
    gene_ids_size = NA,
    geneInfoListLoader = NA
  )
)