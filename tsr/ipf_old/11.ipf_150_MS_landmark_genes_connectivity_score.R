source("modules/config.R")
source("modules/TSRLogger.R")
source("modules/connectivity_score/ConnectivityScoreByGeneSelectionBinChen.R")
source("modules/connectivity_score/OverallConnectivityScore.R")

# setup
tsrLogger <- TSRLogger$new()
overallConnectivityScore <- OverallConnectivityScore$new(ConnectivityScoreByGeneSelectionBinChen$new())
tsrLogger$log("start 150 most significant genes connectivity score computation")

cs_parameters <- list()
cs_parameters$disease_name <- "ipf"
cs_parameters$disease_signature_filename <- gsub(".csv", ".Rds", ipf_disease_signature_cfg$ipf_150_most_significant_landmark_genes_filename)
cs_parameters$drug_signatures_filename <- paste0(drug_signature_cfg$signature_base_dir, ipf_disease_signature_cfg$metanalysis_150_MS_landmark_base_filename, ".Rds")
cs_parameters$drugs_filter <- NA
cs_parameters$perturbation_time_list <- c("6h", "24h", "6h_24h")
cs_parameters$gene_selection_strategy <- disease_signature_cfg$gene_select_strategy_150_most_significant
cs_parameters$n_permutations <- 10^5

connectivity_score <- overallConnectivityScore$compute(cs_parameters)

filename_no_ext <- paste0(connectivity_score_cfg$cs_base_dir, ipf_disease_signature_cfg$connectivity_score_150_MS_landmark_filename)
saveRDS(connectivity_score, paste0(filename_no_ext, ".Rds"))
write.table(connectivity_score, paste0(filename_no_ext, ".csv"), sep = ";", row.names = F, dec = ",")
tsrLogger$log("end 150 most significant genes connectivity score computation")