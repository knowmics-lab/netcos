library(R6)
source("modules/config.R")
source("modules/TSRLogger.R")
source("modules/obj_is_na.R")
source("modules/connectivity_score/merger/Metanalysis6h24hMerger.R")

DrugSignatureLoader <- R6Class(
  "DrugSignatureLoader",
  public = list(
    load = function(disease_genes, drugs_filter = NA) {
      start_time <- Sys.time()
      private$tsrLogger$log("start drug signatures loading")
      drug_signatures <- private$metanalysis6h24hMerger$merge(private$metanalysis_signature_base_dir, disease_genes)
      if (!obj_is_na(drugs_filter)) {
        drug_signatures <- subset(drug_signatures, drug %in% drugs_filter)
      }
      total_time <- Sys.time() - start_time
      private$tsrLogger$log(sprintf("end drug signatures loading, time: %s %s", total_time, attr(total_time, "units")))
      return(drug_signatures)
    }
  ),
  private = list(
    tsrLogger = TSRLogger$new(),
    metanalysis_signature_base_dir = drug_signature_cfg$metanalysis_6h_24h_full_path_dir,
    metanalysis6h24hMerger = Metanalysis6h24hMerger$new()
  )
)