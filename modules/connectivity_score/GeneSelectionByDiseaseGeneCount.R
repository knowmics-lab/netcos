library(R6)
GeneSelectionByDiseaseGeneCount <- R6Class(
  "GeneSelectionByDiseaseGeneCount",
  public = list(
    compute = function(disease_gene_count) {
      gene_selection <- data.frame(
        name = c("50_most_significant", "100_most_significant", "150_most_significant"),
        count = c(50, 100, 150)
      )
      gene_selection <- subset(gene_selection, disease_gene_count + 50 > count)
      if (disease_gene_count < 150) {
        gene_selection[nrow(gene_selection), "count"] <- disease_gene_count
      }
      return(gene_selection)
    }
  )
)