library(R6)
source("modules/config.R")
GeneFilterByProteinCoding <- R6Class(
  "GeneFilterByProteinCoding",
  public = list(
    filter = function(gene_expressions) {
      protein_coding_genes <- readRDS(gene_meta_info_cfg$filename)
      protein_coding_genes <- subset(protein_coding_genes, protein_coding == T)
      gene_expressions$gene <- rownames(gene_expressions)
      gene_expressions <- subset(gene_expressions, gene %in% protein_coding_genes$"Gene name")
      gene_expressions$gene <- NULL
      return(gene_expressions)
    }
  )
)

