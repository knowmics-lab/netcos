library(R6)

source("modules/drug_signature/loader/ExperimentDataLoader.R")
RDS_LINCSExperimentDataLoader <- R6Class(
  "RDS_LINCSExperimentDataLoader",
  inherit = ExperimentDataLoader,
  public = list(
    initialize = function(lincs_splitted_level3_dir) {
      private$lincs_splitted_level3_dir <- lincs_splitted_level3_dir
    },
    load = function(gene_id, experiments_env, computation_number) {
      # legge l'epsressione genica del gene "gene" ottenuta mediante gli esperimenti
      # experiment_info_list$inst_id (sono tanti esperimenti, nell'ordine delle decina di migliaia).
      # Si ottiene una matrice di una riga. La matrice ha come intestazione delle colonne
      # gli id degli esperimenti e come intestazione dell'unica riga l'id del gene. Esempio:
      #       ASG001_MCF7_6H_X2_B7_DUO52HI53LO:A13  ASG001_MCF7_6H_X2_B7_DUO52HI53LO:D13  ASG001_MCF7_6H_X2_B7_DUO52HI53LO:F13
      # 7849  4.581                                 4.0124                                4.2052

      filename <- paste0(private$lincs_splitted_level3_dir, gene_id, ".Rds")
      gene_expressions <- readRDS(filename)
      gene_expressions <- gene_expressions[, experiments_env$experiments_meta_data$inst_id]
      return(gene_expressions)
    }
  ),
  private = list(
    lincs_splitted_level3_dir = NA
  )
)
