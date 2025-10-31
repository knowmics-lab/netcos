library(R6)

source("modules/config.R")
source("modules/obj_is_na.R")
source("modules/TSRLogger.R")
source("modules/drug_signature/DrugSignatureInputArguments.R")
source("modules/drug_signature/loader/GeneSymbolsLoader.R")
source("modules/drug_signature/factory/OverallDrugSignatureFactory.R")

OverallDrugSignatureMain <- R6Class(
  "OverallDrugSignatureMain",
  public = list(
    initialize = function(default_genes_filename = NA) {
      private$drugSignatureInputArguments <- DrugSignatureInputArguments$new(default_genes_filename)
      private$overallDrugSignatureFactory <- OverallDrugSignatureFactory$new()
      private$geneSymbolsLoader <- GeneSymbolsLoader$new()
      private$tsrLogger <- TSRLogger$new()
    },
    compute = function() {
      private$tsrLogger$log("start processing LINCS data")
      inputArguments <- private$drugSignatureInputArguments$get()
      if (obj_is_na(inputArguments$genes_filename)) {
        stop("missing genes file name")
      }
      overallDrugSignature <- private$overallDrugSignatureFactory$create(inputArguments$BLAS_num_threads)
      private$tsrLogger$log(sprintf("gene list file: %s", inputArguments$genes_filename))
      perturbation_times_label <- paste0(paste(inputArguments$perturbation_times, collapse = "h, "), "h")
      private$tsrLogger$log(sprintf("perturbation times: %s", perturbation_times_label))
      gene_symbols <- private$geneSymbolsLoader$load(inputArguments$genes_filename)
      private$tsrLogger$log(sprintf("total genes to compute: %s", length(gene_symbols) * length(inputArguments$perturbation_times)))
      drugs_filter <- read.table(drug_signature_cfg$drug_list, header = T)$drug
      overallDrugSignature$compute(gene_symbols, drugs_filter, inputArguments$perturbation_times)
      private$tsrLogger$log("end processing LINCS data")
      return(NA)
    }
  ),
  private = list(
    drugSignatureInputArguments = NA,
    overallDrugSignatureFactory = NA,
    geneSymbolsLoader = NA,
    tsrLogger = NA
  )
)