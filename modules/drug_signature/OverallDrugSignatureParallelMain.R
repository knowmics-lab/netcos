library(R6)

source("modules/config.R")
source("modules/obj_is_na.R")
source("modules/TSRLogger.R")
source("modules/drug_signature/factory/OverallDrugSignatureFactory.R")

OverallDrugSignatureParallelMain <- R6Class(
  "OverallDrugSignatureParallelMain",
  public = list(
    compute = function(gene_symbols, perturbation_times, BLAS_num_threads) {
      private$tsrLogger$log("start processing LINCS data")
      overallDrugSignature <- private$overallDrugSignatureFactory$create(BLAS_num_threads)
      perturbation_times_label <- paste0(paste(perturbation_times, collapse = "h, "), "h")
      private$tsrLogger$log(sprintf("perturbation times: %s", perturbation_times_label))
      private$tsrLogger$log(sprintf("total genes to compute: %s", length(gene_symbols) * length(perturbation_times)))
      drugs_filter <- read.table(drug_signature_cfg$drug_list, header = T)$drug
      overallDrugSignature$compute(gene_symbols, drugs_filter, perturbation_times)
      private$tsrLogger$log("end processing LINCS data")
      return(NA)
    }
  ),
  private = list(
    overallDrugSignatureFactory = OverallDrugSignatureFactory$new(),
    tsrLogger = TSRLogger$new()
  )
)