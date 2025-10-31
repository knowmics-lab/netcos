library(R6)
source("modules/TSRLogger.R")
source("modules/drug_signature/loader/GeneSymbolsLoader.R")
source("modules/drug_signature/merger/Metanalysis6h24hByGeneSymbolsMerger.R")

Metanalysis6h24hMerger <- R6Class(
  "Metanalysis6h24hMerger",
  public = list(
    merge = function(gene_symbols_filename, metanalysis_signature_base_dir, metanalysis_filename) {
      gene_symbols <- private$geneSymbolsLoader$load(gene_symbols_filename)
      start_time <- Sys.time()
      private$tsrLogger$log(sprintf("start metanalysis merge, total genes: %s", length(gene_symbols)))
      merged_signatures <- private$
        metanalysis6h24hByGeneSymbolsMerger$
        merge(metanalysis_signature_base_dir, gene_symbols)
      saveRDS(merged_signatures, metanalysis_filename)
      metanalysis_filename_csv <- gsub(".Rds", ".csv", metanalysis_filename,)
      write.table(merged_signatures, metanalysis_filename_csv, sep = ";", row.names = F, dec = ",")
      total_time <- Sys.time() - start_time
      private$tsrLogger$log(sprintf("end metanalysis merge, time: %s %s", total_time, attr(total_time, "units")))
      return(NA)
    }
  ),
  private = list(
    geneSymbolsLoader = GeneSymbolsLoader$new(),
    metanalysis6h24hByGeneSymbolsMerger = Metanalysis6h24hByGeneSymbolsMerger$new(),
    tsrLogger = TSRLogger$new()
  )
)