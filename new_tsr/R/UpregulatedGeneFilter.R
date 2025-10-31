
UpregulatedGeneFilter <- R6Class(
  "UpregulatedGeneFilter",
  public = list(
    filter = function(signature) {
      upregulated_genes <- subset(signature, estimate > 0)
      upregulated_genes$gene_id <- rownames(upregulated_genes)
      upregulated_genes <- upregulated_genes[order(upregulated_genes$estimate, decreasing = T),]
      return(upregulated_genes[, "gene_id", drop = F])
    }
  )
)