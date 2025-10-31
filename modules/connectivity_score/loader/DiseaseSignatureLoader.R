library(R6)
DiseaseSignatureLoader <- R6Class(
  "DiseaseSignatureLoader",
  public = list(
    load = function(disease_signature_filename) {
      disease_signature <- readRDS(disease_signature_filename)
      disease_signature <- disease_signature[order(abs(disease_signature$t.value), decreasing = T),]
    }
  )
)