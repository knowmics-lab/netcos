
RandomConnectivityScoreDistribution <- R6Class(
  "RandomConnectivityScoreDistribution",
  inherit = RandomConnectivityScoreDistributionAbstract,
  public = list(
    initialize = function() {
      private$cMapScore <- CMapScore$new()
    },
    compute = function(n_disease_up_regulated_genes, n_disease_down_regulated_genes, n_drug_signatures_genes, n_permutations) {
      output <- numeric(n_permutations)

      for (i in 1:n_permutations) {

        drug_signature <- data.frame(
          ids = 1:n_drug_signatures_genes,
          rank = sample(1:n_drug_signatures_genes, replace = F)
        )

        DEG_genes <- sample(1:(n_disease_up_regulated_genes + n_disease_down_regulated_genes), replace = F)
        sig_up <- DEG_genes[1:n_disease_up_regulated_genes]
        sig_down <- DEG_genes[(n_disease_up_regulated_genes + 1):length(DEG_genes)]

        output[i] <- private$cMapScore$compute(
          sig_up,
          sig_down,
          drug_signature
        )
      }
      return(output)
    }
  ),
  private = list(
    cMapScore = NA
  )
)
