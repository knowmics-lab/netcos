library(R6)

UpregulatedGeneFilter <- R6Class(
  "UpregulatedGeneFilter",
  public = list(
    filter = function(disease_signature) {
      upregulated_genes <- subset(disease_signature, estimate > 0)
      upregulated_genes$GeneID <- rownames(upregulated_genes)
      upregulated_genes <- upregulated_genes[order(upregulated_genes$estimate, decreasing = T),]
      return(upregulated_genes[, "GeneID", drop = F])
    }
  )
)