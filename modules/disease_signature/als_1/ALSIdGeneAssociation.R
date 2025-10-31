library(R6)
library(GEOquery)

ALSIdGeneAssociation <- R6Class(
  "ALSIdGeneAssociation",
  public = list(
    association = function(gse_3307_gpl97) {
      id_gene_association <- fData(gse_3307_gpl97)
      id_gene_association <- id_gene_association[, c("ID", "Gene symbol")]
      colnames(id_gene_association)[2] <- "gene"
      return(id_gene_association)
    }
  )
)