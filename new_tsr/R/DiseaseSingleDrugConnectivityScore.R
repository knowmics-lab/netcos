
DiseaseSingleDrugConnectivityScore <- R6Class(
  "DiseaseSingleDrugConnectivityScore",
  public = list(
    initialize = function(parallel_computation = F) {
      private$diseaseSignatureEstimateMapper <- DiseaseSignatureEstimateMapper$new()
      private$binChenConnectivityScore <- BinChenConnectivityScore$new()
      private$downregulatedGeneFilter <- DownregulatedGeneFilter$new()
      private$upregulatedGeneFilter <- UpregulatedGeneFilter$new()
      private$drugSignatureMapper <- DrugSignatureMapper$new()
      if (parallel_computation) {
        private$randomConnectivityScoreDistribution <- RandomConnectivityScoreDistributionParallel$new()
      }else {
        private$randomConnectivityScoreDistribution <- RandomConnectivityScoreDistribution$new()
      }
      private$connectivityScorePValue <- ConnectivityScorePValue$new()
    },
    compute = function(disease_dge, drug_dge, n_most_significant_genes, t_value_column_name = "t.value", n_permutations = 10^5, compute_p_value = F) {
      disease_dge <- subset(disease_dge, disease_dge$gene_id %in% drug_dge$gene_id)
      drug_signature <- subset(drug_dge, drug_dge$gene_id %in% disease_dge$gene_id)
      drug_signature <- private$drugSignatureMapper$map(drug_signature, t_value_column_name)
      disease_dge$abs_t.value <- abs(disease_dge$t.value)
      disease_dge <- disease_dge[order(disease_dge$abs_t.value, decreasing = T),]
      disease_signature <- disease_dge[1:n_most_significant_genes, , drop = FALSE]
      disease_signature <- private$diseaseSignatureEstimateMapper$map(disease_signature)
      disease_up_regulated_genes <- private$upregulatedGeneFilter$filter(disease_signature)
      disease_down_regulated_genes <- private$downregulatedGeneFilter$filter(disease_signature)
      connectivity_score <- private$
        binChenConnectivityScore$
        compute(disease_down_regulated_genes, disease_up_regulated_genes, drug_signature)
      if (compute_p_value) {
        if ("RandomConnectivityScoreDistributionParallel" %in% class(private$randomConnectivityScoreDistribution)) {
          processorCores$initCores()
        }
        random_connectivity_score_distribution <-
          private$
            randomConnectivityScoreDistribution$
            compute(
            n_disease_up_regulated_genes = dim(disease_up_regulated_genes)[1],
            n_disease_down_regulated_genes = dim(disease_down_regulated_genes)[1],
            n_drug_signatures_genes = dim(drug_signature)[1],
            n_permutations = n_permutations # 10^5
          )
        p.value <- private$connectivityScorePValue$compute(connectivity_score, random_connectivity_score_distribution)
        connectivity_score <- c(connectivity_score, p.value)
      }
      return(connectivity_score)
    }
  ),
  private = list(
    diseaseSignatureEstimateMapper = NA,
    downregulatedGeneFilter = NA,
    upregulatedGeneFilter = NA,
    drugSignatureMapper = NA,
    binChenConnectivityScore = NA,
    randomConnectivityScoreDistribution = NA,
    connectivityScorePValue = NA
  )
)
