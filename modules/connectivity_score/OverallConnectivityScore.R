library(R6)

source("modules/TSRLogger.R")
source("modules/connectivity_score/loader/SignaturesLoader.R")
source("modules/connectivity_score/ConnectivityScoreByPerturbationTime.R")

OverallConnectivityScore <- R6Class(
  "OverallConnectivityScore",
  public = list(
    initialize = function(connectivityScoreByGeneSelection) {
      private$tsrLogger <- TSRLogger$new()
      private$signaturesLoader <- SignaturesLoader$new()
      private$connectivityScoreByPerturbationTime <- ConnectivityScoreByPerturbationTime$new(connectivityScoreByGeneSelection)
    },
    compute = function(cs_parameters) {
      connectivity_score <- data.frame()
      startTime <- Sys.time()
      private$tsrLogger$log(paste0("start ", cs_parameters$computation_title))
      signatures <- private$signaturesLoader$load(cs_parameters$disease_most_significant_genes_filename, cs_parameters$drugs_filter)
      for (perturbation_time in cs_parameters$perturbation_time_list) {
        tmp <- private$
          connectivityScoreByPerturbationTime$
          compute(
          cs_parameters$disease_name,
          signatures$disease,
          signatures$drugs,
          perturbation_time,
          cs_parameters$gene_selection_strategy,
          cs_parameters$n_permutations
        )
        connectivity_score <- rbind(connectivity_score, tmp)
      }
      totalTime <- Sys.time() - startTime
      private$tsrLogger$log(sprintf("end %s, time: %s %s", cs_parameters$computation_title, totalTime, attr(totalTime, "units")))

      connectivity_score_rds_filename <- gsub(".csv", ".Rds", cs_parameters$connectivity_score_filename)
      saveRDS(connectivity_score, connectivity_score_rds_filename)
      write.table(connectivity_score, cs_parameters$connectivity_score_filename, sep = ";", row.names = F, dec = ",")
      return(connectivity_score)
    }
  ),
  private = list(
    tsrLogger = NA,
    signaturesLoader = NA,
    connectivityScoreByPerturbationTime = NA
  )
)
