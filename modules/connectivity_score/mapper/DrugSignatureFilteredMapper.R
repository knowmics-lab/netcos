library(R6)

DrugSignatureFilteredMapper <- R6Class(
  "DrugSignatureFilteredMapper",
  public = list(
    map = function(drug_signatures, drug_name_filter) {
      drug_signatures <- subset(drug_signatures, drug == drug_name_filter)
      drug_signatures <- drug_signatures[!duplicated(drug_signatures$gene),]
      rownames(drug_signatures) <- drug_signatures$gene
      drug_signatures <- drug_signatures[, "t.value", drop = F]
      colnames(drug_signatures) <- "estimate"
      return(drug_signatures)
    }
  )
)