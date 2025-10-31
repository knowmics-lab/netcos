library(R6)
library(cmapR)

source("modules/drug_signature/loader/ExperimentDataLoader.R")
GCTX_LINCSExperimentDataLoader <- R6Class(
  "GCTX_LINCSExperimentDataLoader",
  inherit = ExperimentDataLoader,
  public = list(
    initialize = function(experiments_data_filename) {
      private$experiments_data_filename <- experiments_data_filename
    },

    load = function(gene_id, experiments_env, computation_number) {
      startTime <- Sys.time()
      # 65 GB file, download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
      # legge l'epsressione genica del gene "gene" ottenuta mediante gli esperimenti
      # experiment_info_list$inst_id (sono tanti esperimenti, nell'ordine delle decina di migliaia).
      # Si ottiene una matrice di una riga. La matrice ha come intestazione delle colonne
      # gli id degli esperimenti e come intestazione dell'unica riga l'id del gene. Esempio:
      #       ASG001_MCF7_6H_X2_B7_DUO52HI53LO:A13  ASG001_MCF7_6H_X2_B7_DUO52HI53LO:D13  ASG001_MCF7_6H_X2_B7_DUO52HI53LO:F13
      # 7849  4.581                                 4.0124                                4.2052

      # usa l'unica riga ottenuta dalla chiamata parse_gctx
      gene_expressions <-
        parse_gctx(
          fname = private$experiments_data_filename,
          rid = gene_id,
          cid = experiments_env$experiments_meta_data$inst_id
        )@mat[1,]
      totalTime <- Sys.time() - startTime
      print(sprintf("computation number: %s, end parsing gctx: %s %s ", computation_number, totalTime, attr(totalTime, "units")))
      return(gene_expressions)
    }
  ),
  private = list(
    experiments_data_filename = NA
  )
)
