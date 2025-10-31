
DrugSignatureLoaderAbstract <- R6Class(
  "DrugSignatureLoaderAbstract",
  public = list(
    load = function(drug_name) {
      stop("I'm an abstract method, implement me")
    },
    init = function(disease_drug_common_genes) {
      private$disease_drug_common_genes <- disease_drug_common_genes
    }
  ),
  private = list(
    disease_drug_common_genes = NA
  )
)