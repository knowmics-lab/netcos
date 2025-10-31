
DiseaseDrugConnectivityScoreByBinChenParallel <- R6Class(
  "DiseaseDrugConnectivityScoreByBinChenParallel",
  inherit = DiseaseDrugConnectivityScoreAbstract,
  public = list(
    initialize = function(drugSignatureLoader = NA) {
      private$connectivityScore <- DiseaseDrugConnectivityScoreByBinChen$new(RandomConnectivityScoreDistributionParallel$new(), DrugListConnectivityScoreByBinChenParallel$new(drugSignatureLoader))
    },
    compute = function(disease_signature, drugs, n_drug_signatures_genes, n_permutations = 10^5, disease_name = NA, gene_selection_strategy = NA, drug_perturbation_time = NA) {
      # disease_signature example:
      #         estimate
      # 780     -13,1465
      # 2197    -12,9753
      # 51493   10,4324
      # estimate = t value, row names = gene id

      # drugs example:
      # name         filename
      # ethisterone  drugs/ethisterone.Rds
      # ethoprop     drugs/ethoprop.Rds
      # ethotoin     drugs/ethotoin.Rds

      return(private$connectivityScore$compute(disease_signature, drugs, n_drug_signatures_genes, n_permutations, disease_name, gene_selection_strategy, drug_perturbation_time))
    }
  ),
  private = list(
    connectivityScore = NA
  )
)