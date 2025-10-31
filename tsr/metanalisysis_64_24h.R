source("als_2/config/als_2_signature_cfg.R")
source("modules/TSRLogger.R")
source("modules/drug_signature/OverallMetanalysis6h24hDrugSignature.R")

# setup
tsrLogger <- TSRLogger$new()
overallMetanalysisDrugSignature6h24h <- OverallMetanalysis6h24hDrugSignature$new()

# computation
tsrLogger$log("start all genes meta analysis")
no_print_output <- overallMetanalysisDrugSignature6h24h$compute("data/LINCS-GSE92742/bing_gene_symbols.csv")
tsrLogger$log("end all genes meta analysis")