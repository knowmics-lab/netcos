
DiseaseDrugConnectivityScoreAbstract <- R6Class(
  "DiseaseDrugConnectivityScoreAbstract",
  public = list(
    compute = function(disease_signature, drugs, n_drug_signatures_genes, n_permutations, disease_name, gene_selection_strategy, drug_perturbation_time) {
      stop("I'm an abstract method, implement me")
    }
  )
)