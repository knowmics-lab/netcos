library(R6)
library(metafor)

source("modules/TSRLogger.R")

Metanalysis6h24hDrugSignature <- R6Class(
  "Metanalysis6h24hDrugSignature",
  public = list(
    compute = function(drug_by_gene_signatures) {
      ma_output <- rma(yi = c(drug_by_gene_signatures$DE_log2_FC_6h, drug_by_gene_signatures$DE_log2_FC_24h), sei = c(drug_by_gene_signatures$std.error_6h, drug_by_gene_signatures$std.error_24h))
      drug_by_gene_signatures[c("DE_log2_FC_6h_24h", "std.error_6h_24h", "t.value_6h_24h", "p.value_6h_24h")] <-
        c(ma_output$b[, 1], ma_output$se, ma_output$zval, ma_output$pval)
      return(drug_by_gene_signatures)
    }
  ),
  private = list(
    tsrLogger = TSRLogger$new()
  )
)