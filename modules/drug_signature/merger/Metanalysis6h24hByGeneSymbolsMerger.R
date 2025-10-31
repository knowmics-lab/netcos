library(R6)
source("modules/config.R")
source("modules/drug_signature/loader/GeneInfoListLoader.R")
source("modules/drug_signature/converter/GeneSymbolToIdConverter.R")

Metanalysis6h24hByGeneSymbolsMerger <- R6Class(
  "Metanalysis6h24hByGeneSymbolsMerger",
  public = list(
    initialize = function() {
      geneInfoListLoader <- GeneInfoListLoader$new(drug_signature_cfg$gene_info_filename)
      private$geneSymbolToIdConverter <- GeneSymbolToIdConverter$new(geneInfoListLoader)
    },
    merge = function(metanalysis_signature_base_dir, gene_symbols) {
      gene_ids <- private$geneSymbolToIdConverter$convert(gene_symbols)
      merged_signatures <- data.frame()
      for (i in 1:length(gene_symbols)) {
        signature_filename <- paste0(metanalysis_signature_base_dir, gene_symbols[i], "_", gene_ids[i], "_metanalysis.Rds")
        signature <- readRDS(signature_filename)
        merged_signatures <- rbind(merged_signatures, signature)
      }
      return(merged_signatures)
    }
  ),
  private = list(
    geneSymbolToIdConverter = NA
  )
)