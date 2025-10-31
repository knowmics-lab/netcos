source("als_2/config/als_2_signature_cfg.R")
source("modules/TSRLogger.R")
source("modules/drug_signature/OverallMetanalysis6h24hDrugSignature.R")

# setup
tsrLogger <- TSRLogger$new()
overallMetanalysisDrugSignature6h24h <- OverallMetanalysis6h24hDrugSignature$new()

# computation
tsrLogger$log("start 150 most significant genes meta analysis")
no_print_output <- overallMetanalysisDrugSignature6h24h$compute(als_2_signature_cfg$a150_most_significant_genes_filename)
tsrLogger$log("end 150 most significant genes meta analysis")
tsrLogger$log("start 150 most significant landmark genes meta analysis")
no_print_output <- overallMetanalysisDrugSignature6h24h$compute(als_2_signature_cfg$a150_most_significant_landmark_genes_filename)
tsrLogger$log("end 150 most significant landmark genes meta analysis")