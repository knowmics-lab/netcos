library(R6)

source("modules/TSRLogger.R")
source("modules/connectivity_score/GeneSelectionStrategy.R")
source("modules/connectivity_score/filter/DrugSignatureFilter.R")

ConnectivityScoreByPerturbationTime <- R6Class(
  "ConnectivityScoreByPerturbationTime",
  public = list(
    initialize = function(connectivityScoreByGeneSelection) {
      if (!"ConnectivityScoreByGeneSelection" %in% class(connectivityScoreByGeneSelection))
        stop("the connectivityScoreByGeneSelection class type must be a subclass of ConnectivityScoreByGeneSelection")
      private$tsrLogger <- TSRLogger$new()
      private$geneSelectionStrategy <- GeneSelectionStrategy$new()
      private$drugSignatureFilter <- DrugSignatureFilter$new()
      private$connectivityScoreByGeneSelection <- connectivityScoreByGeneSelection
    },
    compute = function(disease_name, disease_signature, drug_signatures, perturbation_time, gene_selection_strategy, n_permutations = 10^5) {
      startTime <- Sys.time()
      private$tsrLogger$log(sprintf("start connectivity score computation by perturbation time %s", perturbation_time))
      # carica tutte le firme della malattia
      total_disease_genes <- nrow(disease_signature)
      if (total_disease_genes < 10) {
        return(NA)
      }
      gene_selection_list <- private$geneSelectionStrategy$compute(total_disease_genes, gene_selection_strategy)
      connectivity_score <- data.frame()
      total_selection_lists <- nrow(gene_selection_list)
      drug_signatures <- private$drugSignatureFilter$filter(drug_signatures, perturbation_time)
      private$tsrLogger$log(sprintf("total selection lists: %s", total_selection_lists))
      for (i in 1:total_selection_lists) {
        gene_selection <- gene_selection_list[i,]
        disease_signature_selected_genes <- disease_signature[1:(gene_selection$count), , drop = FALSE]
        tmp <- private$connectivityScoreByGeneSelection$compute(disease_name, disease_signature_selected_genes, drug_signatures, gene_selection$name, perturbation_time, n_permutations)
        connectivity_score <- rbind(connectivity_score, tmp)
      }
      totalTime <- Sys.time() - startTime
      private$tsrLogger$log(sprintf("end connectivity score computation by perturbation time %s, time: %s %s", perturbation_time, totalTime, attr(totalTime, "units")))
      return(connectivity_score)
    }
  ),
  private = list(
    tsrLogger = NA,
    geneSelectionStrategy = NA,
    drugSignatureFilter = NA,
    connectivityScoreByGeneSelection = NA
  )
)