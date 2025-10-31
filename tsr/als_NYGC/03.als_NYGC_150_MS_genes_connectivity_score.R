source("als_NYGC/config/als_NYGC_cfg.R")
source("modules/connectivity_score/ConnectivityScoreByGeneSelectionBinChen.R")
source("modules/connectivity_score/OverallConnectivityScore.R")

overallConnectivityScore <- OverallConnectivityScore$new(ConnectivityScoreByGeneSelectionBinChen$new())

cs_parameters <- list()
cs_parameters$computation_title <- "als NYGC 150 most significant genes connectivity score computation"
cs_parameters$disease_name <- "als"
cs_parameters$disease_most_significant_genes_filename <- gsub(".csv", ".Rds", als_NYGC_cfg$signature_150_most_significant_genes_filename)
cs_parameters$drugs_filter <- NA
cs_parameters$perturbation_time_list <- c("6h", "24h", "6h_24h")
cs_parameters$gene_selection_strategy <- disease_signature_cfg$gene_select_strategy_150_most_significant
cs_parameters$n_permutations <- 10^5
cs_parameters$connectivity_score_filename <- als_NYGC_cfg$connectivity_score_150_MS_genes_filename

connectivityScore <- overallConnectivityScore$compute(cs_parameters)
