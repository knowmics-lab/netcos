library(R6)

ConnectivityScoreByGeneSelection <- R6Class(
  "ConnectivityScoreByGeneSelection",
  public = list(
    compute = function(disease_name, disease_signature, drug_signatures, gene_selection, perturbation_time, n_permutations = 10^5) {
      stop("I'm an abstract method, implement me")
    }
  )
)