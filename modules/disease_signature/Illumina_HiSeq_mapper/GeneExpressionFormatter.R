library(R6)
source("modules/disease_signature/Illumina_HiSeq_mapper/IdGeneAssociation.R")
source("modules/disease_signature/AddGeneColumn.R")

GeneExpressionFormatter <- R6Class(
  "GeneExpressionFormatter",
  public = list(
    format = function(gene_expressions) {
      gene_expressions <- data.frame(gene_expressions)
      gene_expressions$ID <- rownames(gene_expressions)
      id_gene_association <- private$idGeneAssociation$load()
      gene_expressions <- private$addGeneColumn$add(gene_expressions, id_gene_association)
      rownames(gene_expressions) <- gene_expressions$gene
      gene_expressions$gene <- NULL
      return(gene_expressions)
    }
  ),
  private = list(
    idGeneAssociation = IdGeneAssociation$new(),
    addGeneColumn = AddGeneColumn$new()
  )
)