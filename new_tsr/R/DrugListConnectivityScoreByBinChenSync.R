
DrugListConnectivityScoreByBinChenSync <- R6Class(
  "DrugListConnectivityScoreByBinChenSync",
  public = list(
    initialize = function(drugSignatureLoader = NA) {
      if (!obj_is_na(drugSignatureLoader)) {
        if (!"DrugSignatureLoaderAbstract" %in% class(drugSignatureLoader))
          stop("the drugSignatureLoader instance must by of type DrugSignatureLoaderAbstract")
        private$drugSignatureLoader <- drugSignatureLoader
      }else {
        private$drugSignatureLoader <- DrugSignatureLoaderByFilename$new()
      }
      private$connectivityScoreByBinChen <- BinChenConnectivityScoreWithPValue$new()
    },
    compute = function(disease_down_regulated_genes, disease_up_regulated_genes, drugs, random_connectivity_score) {
      total_drugs <- dim(drugs)[1]
      connectivity_score_matrix <- data.frame()
      for (i in 1:total_drugs) {
        drug_signature <- private$drugSignatureLoader$load(drugs$filename[i])
        connectivity_score <- private$connectivityScoreByBinChen$compute(disease_down_regulated_genes, disease_up_regulated_genes, drug_signature, random_connectivity_score)
        connectivity_score_matrix <- rbind(connectivity_score_matrix, data.frame(drugs$name[i], connectivity_score[1], connectivity_score[2]))
      }
      return(connectivity_score_matrix)
    }
  ),
  private = list(
    drugSignatureLoader = NA,
    connectivityScoreByBinChen = NA
  )
)