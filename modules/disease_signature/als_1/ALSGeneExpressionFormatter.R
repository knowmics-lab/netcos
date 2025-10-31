library(R6)
source("modules/disease_signature/AddGeneColumn.R")
source("modules/disease_signature/als_1/GeneAggregator.R")
source("modules/disease_signature/als_1/ALSIdGeneAssociation.R")
ALSGeneExpressionFormatter <- R6Class(
  "ALSGeneExpressionFormatter",
  public = list(
    format = function(gene_expressions, gse_3307_gpl97) {
      gene_expressions <- data.frame(gene_expressions)
      id_gene_association <- private$alsIdGeneAssociation$association(gse_3307_gpl97)
      gene_expressions <- private$addGeneColumn$add(gene_expressions, id_gene_association)
      gene_expressions <- private$geneAggregator$aggregate_by_gene(gene_expressions)
    }
  ),
  private = list(
    addGeneColumn = AddGeneColumn$new(),
    geneAggregator = GeneAggregator$new(),
    alsIdGeneAssociation = ALSIdGeneAssociation$new()
  )
)