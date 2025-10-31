source("modules/config.R")
source("modules/connectivity_score/OverallConnectivityScore.R")
source("modules/TSRLogger.R")

# setup
tsrLogger <- TSRLogger$new()
overallConnectivityScore <- OverallConnectivityScore$new()
tsrLogger$log("start connectivity score computation")

cs_parameters <- list()
cs_parameters$disease_name <- "als"
cs_parameters$disease_signature_filename <- paste0(als_1_disease_signature_cfg$signature_base_dir, als_1_disease_signature_cfg$signature_filename, ".Rds")
cs_parameters$drug_signatures_filename <- paste0(drug_signature_cfg$signature_base_dir, drug_signature_cfg$signature_ma_filename, ".Rds")
cs_parameters$drugs_filter <- NA
cs_parameters$perturbation_time_list <- c("6h", "24h", "6h_24h")
cs_parameters$n_permutations <- 10

connectivity_score <- overallConnectivityScore$compute(cs_parameters)

filename_no_ext <- paste0(connectivity_score_cfg$cs_base_dir, connectivity_score_cfg$cs_filename)
saveRDS(connectivity_score, paste0(filename_no_ext, ".Rds"))
write.table(connectivity_score, paste0(filename_no_ext, ".csv"), sep = ";", row.names = F)
print("end connectivity score computation")