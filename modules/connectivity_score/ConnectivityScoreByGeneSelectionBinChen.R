library(R6)
library(foreach)

source("modules/connectivity_score/abstract_class/ConnectivityScoreByGeneSelection.R")
source("modules/TSRLogger.R")
source("modules/CoreProcessorInit.R")
source("modules/connectivity_score/cmap_score_new/RandomConnectivityScoreParallel.R")
source("modules/connectivity_score/filter/DownregulatedGeneFilter.R")
source("modules/connectivity_score/mapper/DrugSignatureFilteredMapper.R")
source("modules/connectivity_score/filter/UpregulatedGeneFilter.R")
source("modules/connectivity_score/cmap_score_new/BinChenConnectivityScore.R")
source("modules/connectivity_score/builder/ConnectivityScoreBuilder.R")

ConnectivityScoreByGeneSelectionBinChen <- R6Class(
  "ConnectivityScoreByGeneSelectionBinChen",
  inherit = ConnectivityScoreByGeneSelection,
  public = list(
    compute = function(disease_name, disease_signature, drug_signatures, gene_selection, perturbation_time, n_permutations = 10^5) {
      # esempio di disease_signature
      #        estimate
      # DDR1   -13,1465
      # FAU    -12,9753
      # RTCB   10,4324

      # genera la distribuzione random necessaria per il cacloclo del p-value
      startTime <- Sys.time()
      private$tsrLogger$log(sprintf("start connectivity score computation by gene selection: %s", gene_selection))
      private$coreProcessorInit$init()
      random_connectivity_score <-
        private$randomConnectivityScoreParallel$compute(
          n_genes_up = sum(disease_signature$estimate > 0),
          n_genes_down = sum(disease_signature$estimate < 0),
          n_genes = length(unique(drug_signatures$gene)),
          n_permutations = n_permutations # 10^5
        )

      drugs <- unique(drug_signatures$drug)
      total_drugs <- length(drugs)
      down_regulated_genes <- private$downregulatedGeneFilter$filter(disease_signature)
      up_regulated_genes <- private$upregulatedGeneFilter$filter(disease_signature)
      connectivity_score_matrix <- foreach(i = 1:total_drugs, .combine = rbind) %dopar% {
        selected_drug_signature <- private$drugSignatureFilteredMapper$map(drug_signatures, drugs[i])
        if (nrow(selected_drug_signature) > 0) {
          connectivity_score_by_drug <- private$binChenConnectivityScore$compute(down_regulated_genes, up_regulated_genes, selected_drug_signature, random_connectivity_score)
          c(i, connectivity_score_by_drug[1], connectivity_score_by_drug[2])
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
    coreProcessorInit = CoreProcessorInit$new(),
    randomConnectivityScoreParallel = RandomConnectivityScoreParallel$new(),
    downregulatedGeneFilter = DownregulatedGeneFilter$new(),
    upregulatedGeneFilter = UpregulatedGeneFilter$new(),
    drugSignatureFilteredMapper = DrugSignatureFilteredMapper$new(),
    binChenConnectivityScore = BinChenConnectivityScore$new(),
    connectivityScoreBuilder = ConnectivityScoreBuilder$new()
  )
)