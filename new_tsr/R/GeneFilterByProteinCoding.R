
GeneFilterByProteinCoding <- R6Class(
  "GeneFilterByProteinCoding",
  public = list(
    initialize = function() {
      private$protein_coding_gene <- package_readRDS(config$protein_coding_gene_filename)
      private$protein_coding_gene <- private$protein_coding_gene[private$protein_coding_gene$locusGroup=='protein-coding gene',]
      private$protein_coding_gene$proteinCoding <- NULL
    },
    filterById = function(gene_expressions) {
      return(gene_expressions[rownames(gene_expressions) %in% private$protein_coding_gene$id,])
    },
    filterBySymbol = function(gene_expressions) {
      return(gene_expressions[rownames(gene_expressions) %in% private$protein_coding_gene$symbol,])
    },
    filterByEnsemblId = function(gene_expressions) {
      return(gene_expressions[rownames(gene_expressions) %in% private$protein_coding_gene$ensemblId,])
    }
  ),
  private = list(
    protein_coding_gene = NA
  )
)

