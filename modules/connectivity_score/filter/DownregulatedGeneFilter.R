library(R6)

DownregulatedGeneFilter <- R6Class(
  "DownregulatedGeneFilter",
  public = list(
    filter = function(signature) {
      downregulated_genes <- subset(signature, estimate < 0)
      downregulated_genes$GeneID <- rownames(downregulated_genes)
      downregulated_genes <- downregulated_genes[order(downregulated_genes$estimate, decreasing = F),]
      return(downregulated_genes[, "GeneID", drop = F])
    }
  )
)