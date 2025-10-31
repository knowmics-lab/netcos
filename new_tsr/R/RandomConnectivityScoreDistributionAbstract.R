
RandomConnectivityScoreDistributionAbstract <- R6Class(
  "RandomConnectivityScoreDistributionAbstract",
  public = list(
    compute = function(n_disease_up_regulated_genes, n_disease_down_regulated_genes, n_drug_signatures_genes, n_permutations) {
      stop("I'm an abstract method, implement me")
    }
  )
)
