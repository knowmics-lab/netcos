library(R6)
source("modules/connectivity_score/loader/DiseaseSignatureLoader.R")
source("modules/connectivity_score/loader/DrugSignatureLoader.R")
source("modules/connectivity_score/mapper/DiseaseSignatureMapper.R")

SignaturesLoader <- R6Class(
  "SignaturesLoader",
  public = list(
    load = function(disease_most_significant_genes_filename, drugs_filter = NA) {
      disease_signature <- private$diseaseSignatureLoader$load(disease_most_significant_genes_filename)
      disease_genes <- unique(disease_signature$gene)
      drug_signatures <- private$drugSignatureLoader$load(disease_genes, drugs_filter)
      signatures <- list()
      signatures$disease <- private$diseaseSignatureMapper$map(disease_signature)
      signatures$drugs <- drug_signatures
      return(signatures)
    }
  ),
  private = list(
    diseaseSignatureLoader = DiseaseSignatureLoader$new(),
    drugSignatureLoader = DrugSignatureLoader$new(),
    diseaseSignatureMapper = DiseaseSignatureMapper$new()
  )
)