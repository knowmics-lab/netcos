library(R6)
library(foreach)

source("modules/drug_signature/Metanalysis6h24hDrugSignatureByGene.R")
source("modules/CoreProcessorInit.R")
source("modules/TSRLogger.R")
source("modules/obj_is_na.R")

OverallMetanalysis6h24hDrugSignature <- R6Class(
  "OverallMetanalysis6h24hDrugSignature",
  public = list(
    initialize = function(drug_by_gene_signates_dir = NA, metanalysis_6h_24_dir = NA) {
      private$metanalysis6h24hDrugSignatureByGene <- Metanalysis6h24hDrugSignatureByGene$new(drug_by_gene_signates_dir, metanalysis_6h_24_dir)
      private$coreProcessorInit <- CoreProcessorInit$new()
      private$tsrLogger <- TSRLogger$new()
    },
    compute = function(gene_list_filename) {
      startTime <- Sys.time()
      private$tsrLogger$log("start meta analysis drug by gene signature 6h/24h")
      gene_list <- read.table(gene_list_filename, sep = ";", header = T)
      total_genes <- dim(gene_list)[1]
      private$tsrLogger$log(sprintf("total genes to compute: %s", total_genes))
      filenames <- paste0(gene_list$gene, "_", gene_list$gene_id)

      private$coreProcessorInit$init()
      foreach(i = 1:total_genes) %dopar% {
        private$
          metanalysis6h24hDrugSignatureByGene$
          compute(filenames[i])
      }
      totalTime <- Sys.time() - startTime
      private$tsrLogger$log(sprintf("end meta analysis drug by gene signature 6h/24h, computation time: %s %s", totalTime, attr(totalTime, "units")))
      return(NA)
    }
  ),
  private = list(
    metanalysis6h24hDrugSignatureByGene = NA,
    coreProcessorInit = NA,
    tsrLogger = NA
  )
)