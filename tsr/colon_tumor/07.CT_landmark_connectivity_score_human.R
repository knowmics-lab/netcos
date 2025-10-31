source("modules/config.R")
source("modules/TSRLogger.R")
source("modules/connectivity_score/OverallConnectivityScore.R")
source("modules/connectivity_score/ConnectivityScoreByGeneSelectionBinChen.R")

# setup
tsrLogger <- TSRLogger$new()
overallConnectivityScore <- OverallConnectivityScore$new(ConnectivityScoreByGeneSelectionBinChen$new())
tsrLogger$log("start colon tumor connectivity score")

cs_parameters <- list()
cs_parameters$disease_name <- "colon tumor human"
cs_parameters$disease_signature_filename <- "data/colon_tumor/colon_tumor_150_most_significant_landmark_genes.Rds"
cs_parameters$drug_signatures_filename <- "output/drug_signature/LINCS/colon_tumor_150_most_significant_landmark_genes_metanalysis_6h_24h.Rds"
cs_parameters$drugs_filter <- NA
cs_parameters$perturbation_time_list <- c("6h", "24h", "6h_24h")
cs_parameters$gene_selection_strategy <- disease_signature_cfg$gene_select_strategy_bin_chen
cs_parameters$n_permutations <- 10^5

connectivity_score <- overallConnectivityScore$compute(cs_parameters)

filename_no_ext <- "output/connectivity_score/colon_tumor_150_most_significant_landmark_genes_connectivity_score"
saveRDS(connectivity_score, paste0(filename_no_ext, ".Rds"))
write.table(connectivity_score, paste0(filename_no_ext, ".csv"), sep = ";", row.names = F, dec = ",")
tsrLogger$log("end colon tumor connectivity score")