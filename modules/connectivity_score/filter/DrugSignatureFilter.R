library(R6)
source("modules/obj_is_na.R")
source("modules/config.R")

DrugSignatureFilter <- R6Class(
  "DrugSignatureFilter",
  public = list(
    filter = function(drug_signatures, perturbation_time) {
      drug_signatures <- drug_signatures[, c("drug", "gene", paste0("t.value_", perturbation_time))]
      colnames(drug_signatures)[3] <- "t.value"
      return(drug_signatures)
    }
  )
)