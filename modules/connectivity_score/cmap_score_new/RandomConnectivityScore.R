library(R6)
source("modules/connectivity_score/cmap_score_new/cmap_score_new.R")

RandomConnectivityScore <- R6Class(
  "RandomConnectivityScore",
  public = list(
    compute = function(n_genes_up, n_genes_down, n_genes, n_permutations) {
      output <- numeric(n_permutations)

      for (i in 1:n_permutations) {

        drug_signature <- data.frame(
          ids = 1:n_genes,
          rank = sample(1:n_genes, replace = F)
        )

        DEG_genes <- sample(1:(n_genes_up + n_genes_down), replace = F)
        sig_up <- DEG_genes[1:n_genes_up]
        sig_down <- DEG_genes[(n_genes_up + 1):length(DEG_genes)]

        output[i] <- cmap_score_new(
          sig_up,
          sig_down,
          drug_signature
        )
      }
      return(output)
    }
  )
)
