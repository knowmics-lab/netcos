library(R6)

source("modules/TSRLogger.R")
source("modules/connectivity_score/abstract_class/ConnectivityScoreByGeneSelection.R")
source("modules/connectivity_score/mapper/DrugSignatureFilteredMapper.R")
source("modules/CoreProcessorInit.R")
source("modules/connectivity_score/builder/ConnectivityScoreBuilder.R")

ConnectivityScoreByGeneSelectionPharmacoGx <- R6Class(
  "ConnectivityScoreByGeneSelectionPharmacoGx",
  inherit = ConnectivityScoreByGeneSelection,
  public = list(
    compute = function(disease_name, disease_signature, drug_signatures, gene_selection, perturbation_time, n_permutations = 10^4) {
      startTime <- Sys.time()
      private$tsrLogger$log(sprintf("start connectivity score computation by gene selection: %s", gene_selection))
      # esempio di disease_signature mappata
      #        estimate
      # DDR1   -13,1465
      # FAU    -12,9753
      # RTCB   10,4324

      drugs <- unique(drug_signatures$drug)
      total_drugs <- length(drugs)
      connectivity_score_matrix <- data.frame();
      for (i in 1:total_drugs) {
        selected_drug_signature <- private$drugSignatureFilteredMapper$map(drug_signatures, drugs[i])
        if (nrow(selected_drug_signature) > 0) {
          connectivity_score_by_drug <- PharmacoGx::connectivityScore(x = selected_drug_signature, y = disease_signature, method = "fgsea", nperm = n_permutations)
          connectivity_score_matrix <- rbind(connectivity_score_matrix, c(i, connectivity_score_by_drug[1], connectivity_score_by_drug[2]))
        }
      }
      connectivity_scores <- private$connectivityScoreBuilder$build(disease_name, drugs, connectivity_score_matrix, gene_selection, perturbation_time)
      totalTime <- Sys.time() - startTime
      private$tsrLogger$log(sprintf("end connectivity score computation by gene selection: %s, time: %s %s", gene_selection, totalTime, attr(totalTime, "units")))
      return(connectivity_scores)
    }
  ),
  private = list(
    tsrLogger = TSRLogger$new(),
    drugSignatureFilteredMapper = DrugSignatureFilteredMapper$new(),
    connectivityScoreBuilder = ConnectivityScoreBuilder$new()
  )
)