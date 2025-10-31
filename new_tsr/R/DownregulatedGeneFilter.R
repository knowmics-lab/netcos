
DownregulatedGeneFilter <- R6Class(
  "DownregulatedGeneFilter",
  public = list(
    filter = function(signature) {
      downregulated_genes <- subset(signature, estimate < 0)
      downregulated_genes$gene_id <- rownames(downregulated_genes)
      downregulated_genes <- downregulated_genes[order(downregulated_genes$estimate, decreasing = F),]
      return(downregulated_genes[, "gene_id", drop = F])
    }
  )
)