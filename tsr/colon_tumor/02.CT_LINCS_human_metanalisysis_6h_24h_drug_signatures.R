source("modules/config.R")
source("modules/drug_signature/OverallMetanalysis6h24hDrugSignature.R")
source("modules/TSRLogger.R")

# setup
tsrLogger <- TSRLogger$new()
overallMetanalysisDrugSignature6h24h <- OverallMetanalysis6h24hDrugSignature$new()

# computation
tsrLogger$log("start human colon tumor genes meta analysis")
filename <- "data/colon_tumor/colon_tumor_150_most_significant_genes_b.csv"
no_print_output <- overallMetanalysisDrugSignature6h24h$compute(filename)
tsrLogger$log("end human colon tumor genes meta analysis")