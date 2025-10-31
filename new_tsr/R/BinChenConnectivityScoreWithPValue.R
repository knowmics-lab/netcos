
BinChenConnectivityScoreWithPValue <- R6Class(
  "BinChenConnectivityScoreWithPValue",
  public = list(
    initialize = function() {
      private$binChenConnectivityScore <- BinChenConnectivityScore$new()
      private$connectivityScorePValue <- ConnectivityScorePValue$new()
    },
    compute = function(disease_down_regulated_genes, disease_up_regulated_genes, drug_signature, random_cs_distribution) {
      connectivity_score <- private$binChenConnectivityScore$compute(disease_down_regulated_genes, disease_up_regulated_genes, drug_signature)
      return(c(connectivity_score, private$connectivityScorePValue$compute(connectivity_score, random_cs_distribution)))
    }
  ),
  private = list(
    binChenConnectivityScore = NA,
    connectivityScorePValue = NA
  )
)






