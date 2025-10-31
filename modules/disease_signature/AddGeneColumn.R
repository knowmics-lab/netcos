library(R6)
AddGeneColumn <- R6Class(
  "AddGeneColumn",
  public = list(
    add = function(gene_expressions, id_gene_association) {
      gene_expressions$ID <- rownames(gene_expressions)
      gene_expressions <- merge(x = id_gene_association, y = gene_expressions, by.x = "ID", by.y = "ID")
      gene_expressions$ID <- NULL
      gene_expressions <- subset(gene_expressions, gene_expressions$gene != "")
      rownames(gene_expressions) <- NULL
      return(gene_expressions)
    }
  )
)