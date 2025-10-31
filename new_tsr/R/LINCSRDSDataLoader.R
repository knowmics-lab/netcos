
LINCSRDSDataLoader <- R6Class(
  "LINCSRDSDataLoader",
  inherit = DataLoaderAbstract,
  public = list(
    initialize = function(lincs_splitted_level3_dir) {
      private$lincs_splitted_level3_dir <- lincs_splitted_level3_dir
    },
    load = function(gene_id, metadata) {
      # legge l'epsressione genica del gene "gene_id" ottenuta mediante gli esperimenti
      # experiment_info_list$inst_id (sono tanti esperimenti, nell'ordine delle decina di migliaia).
      # Si ottiene una matrice di una riga. La matrice ha come intestazione delle colonne
      # gli id degli esperimenti e come intestazione dell'unica riga l'id del gene. Esempio:
      #       ASG001_MCF7_6H_X2_B7_DUO52HI53LO:A13  ASG001_MCF7_6H_X2_B7_DUO52HI53LO:D13  ASG001_MCF7_6H_X2_B7_DUO52HI53LO:F13
      # 7849  4.581                                 4.0124                                4.2052

      filename <- paste0(private$lincs_splitted_level3_dir, gene_id, ".Rds")
      gene_expression <- readRDS(filename)
      metadata$gene_expression <- gene_expression[, metadata$inst_id]
      return(metadata)
    }
  ),
  private = list(
    lincs_splitted_level3_dir = NA
  )
)
