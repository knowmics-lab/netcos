library(R6)

DiseaseSignatureMapper <- R6Class(
  "DiseaseSignatureMapper",
  public = list(
    map = function(disease_signature) {
      rownames(disease_signature) <- disease_signature$gene
      disease_signature <- disease_signature[, "t.value", drop = F]
      colnames(disease_signature) <- "estimate"
      return(disease_signature)
    }
  )
)
