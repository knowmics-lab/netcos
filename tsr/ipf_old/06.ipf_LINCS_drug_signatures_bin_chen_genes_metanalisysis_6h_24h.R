source("modules/config.R")
source("modules/TSRLogger.R")
source("modules/drug_signature/OverallMetanalysis6h24hDrugSignature.R")

# setup
tsrLogger <- TSRLogger$new()
overallMetanalysisDrugSignature6h24h <- OverallMetanalysis6h24hDrugSignature$new()

# computation
tsrLogger$log("start merge ipf bin chen most significant genes drug signature")
no_print_output <- overallMetanalysisDrugSignature6h24h$compute(ipf_disease_signature_cfg$ipf_bin_chen_most_significant_genes_filename)
tsrLogger$log("end merge ipf bin chen most significant genes drug signature")