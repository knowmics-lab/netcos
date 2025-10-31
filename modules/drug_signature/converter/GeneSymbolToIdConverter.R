library(R6)

GeneSymbolToIdConverter <- R6Class(
  "GeneSymbolToIdConverter",
  public = list(
    initialize = function(geneInfoListLoader) {
      if (!"GeneInfoListLoader" %in% class(geneInfoListLoader))
        stop("geneInfoListLoader must by an istance of GeneInfoListLoader")
      private$geneInfoListLoader <- geneInfoListLoader
    },

    convert = function(gene_symbol_list) {
      gene_info_list <- private$geneInfoListLoader$load()
      gene_info_list <- subset(gene_info_list, gene_info_list$pr_gene_symbol %in% gene_symbol_list)
      gene_info_list <- gene_info_list[match(gene_symbol_list,gene_info_list$pr_gene_symbol ),]
      return(as.character(gene_info_list$pr_gene_id))
    }
  ),
  private = list(
    geneInfoListLoader = NA
  )
)