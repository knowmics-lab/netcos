source("modules/config.R")
source("modules/drug_signature/OverallMetanalysis6h24hDrugSignature.R")
source("modules/TSRLogger.R")

# setup
tsrLogger <- TSRLogger$new()
overallMetanalysisDrugSignature6h24h <- OverallMetanalysis6h24hDrugSignature$new()

# computation
tsrLogger$log("start 150 most significant genes meta analysis")
no_print_output <- overallMetanalysisDrugSignature6h24h$compute(ipf_disease_signature_cfg$ipf_150_most_significant_genes_filename)
tsrLogger$log("end 150 most significant genes meta analysis")