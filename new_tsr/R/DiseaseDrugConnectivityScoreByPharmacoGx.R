
DiseaseDrugConnectivityScoreByPharmacoGx <- R6Class(
  "DiseaseDrugConnectivityScoreByPharmacoGx",
  inherit = DiseaseDrugConnectivityScoreAbstract,
  public = list(
    initialize = function(drugSignatureLoader = NA) {
      if (obj_is_na(drugSignatureLoader)) {
        private$drugSignatureLoader <- DrugSignatureLoaderByFilename$new()
      }else {
        private$drugSignatureLoader <- drugSignatureLoader
      }
      private$drugListConnectivityScoreMapper <- DrugListConnectivityScoreMapper$new()
    },
    compute = function(disease_signature, drugs, n_drug_signatures_genes = NA, n_permutations = 10^4, disease_name = NA, gene_selection_strategy = NA, drug_perturbation_time = NA) {
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

      startTime <- Sys.time()
      dgrpLogger$log(sprintf("start connectivity score computation by PharmacoGx Algorithm"))
      total_drugs <- dim(drugs)[1]
      connectivity_score_matrix <- data.frame()
      for (i in 1:total_drugs) {
        drug_signature <- private$drugSignatureLoader$load(drugs$filename[i])
        connectivity_score_by_drug <- PharmacoGx::connectivityScore(x = drug_signature, y = disease_signature, method = "fgsea", nperm = n_permutations)
        connectivity_score_matrix <- rbind(connectivity_score_matrix, data.frame(drugs$name[i], connectivity_score_by_drug[1], connectivity_score_by_drug[2]))
      }
      connectivity_scores <- private$drugListConnectivityScoreMapper$map(connectivity_score_matrix, disease_name, gene_selection_strategy, drug_perturbation_time)
      totalTime <- Sys.time() - startTime
      dgrpLogger$log(sprintf("end connectivity score computation by PharmacoGx Algorithm: time: %s %s", totalTime, attr(totalTime, "units")))
      return(connectivity_scores)
    }
  ),
  private = list(
    drugSignatureLoader = NA,
    drugListConnectivityScoreMapper = NA
  )
)
