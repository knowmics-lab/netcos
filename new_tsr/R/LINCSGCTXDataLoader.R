
LINCSGCTXDataLoader <- R6Class(
  "LINCSGCTXDataLoader",
  inherit = DataLoaderAbstract,
  public = list(
    initialize = function(lincs_level3_filename) {
      private$lincs_level3_filename <- lincs_level3_filename
    },

    load = function(gene_id, metadata) {
      startTime <- Sys.time()
      # 65 GB file, download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
      # legge l'epsressione genica del gene "gene_id" ottenuta mediante gli esperimenti
      # experiment_info_list$inst_id (sono tanti esperimenti, nell'ordine delle decina di migliaia).
      # Si ottiene una matrice di una riga. La matrice ha come intestazione delle colonne
      # gli id degli esperimenti e come intestazione dell'unica riga l'id del gene. Esempio:
      #       ASG001_MCF7_6H_X2_B7_DUO52HI53LO:A13  ASG001_MCF7_6H_X2_B7_DUO52HI53LO:D13  ASG001_MCF7_6H_X2_B7_DUO52HI53LO:F13
      # 7849  4.581                                 4.0124                                4.2052

      # usa l'unica riga ottenuta dalla chiamata parse_gctx
      dgrpLogger$log(sprintf("start reading gctx archive by gene id: %s", gene_id))
      gene_expression <-
        parse_gctx(
          fname = private$lincs_level3_filename,
          rid = gene_id,
          cid = metadata$inst_id
        )@mat[1,]
      totalTime <- Sys.time() - startTime
      dgrpLogger$log(sprintf("end reading gctx archive by gene id: %s, time: %s %s ", gene_id, totalTime, attr(totalTime, "units")))
      metadata$gene_expression <- gene_expression
      return(metadata)
    }
  ),
  private = list(
    lincs_level3_filename = NA
  )
)
