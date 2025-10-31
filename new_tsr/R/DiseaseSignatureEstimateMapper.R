
DiseaseSignatureEstimateMapper <- R6Class(
  "DiseaseSignatureEstimateMapper",
  public = list(
    map = function(disease_signature) {
      rownames(disease_signature) <- disease_signature$gene_id
      disease_signature <- disease_signature[, "t.value", drop = F]
      colnames(disease_signature) <- "estimate"
      return(disease_signature)
    }
  )
)
