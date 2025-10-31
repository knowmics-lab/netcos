source("modules/config.R")
source("modules/drug_signature/OverallDrugSignature.R")
source("modules/drug_signature/OverallMetanalysisDrugSignature6h24h.R")
source("modules/drug_signature/GeneInfoListLoader.R")
source("modules/drug_signature/LINCSExperimentMetaDataSetuper.R")
source("modules/drug_signature/GCTX_LINCSExperimentDataLoader.R")

# setup
geneInfoListLoader <- GeneInfoListLoader$new(drug_signature_cfg$gene_info_filename)
lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)
gctx_LINCSExperimentDataLoader <- GCTX_LINCSExperimentDataLoader$new(drug_signature_cfg$experiments_data_filename)
overallDrugSignature <- OverallDrugSignature$new(geneInfoListLoader, lincsExperimentMetaDataSetuper, gctx_LINCSExperimentDataLoader)
overallMetanalysisDrugSignature6h24h <- OverallMetanalysisDrugSignature6h24h$new()


print("start processing LINCS data")

gene_symbols <- c("DDR1", "PAX8", "FAU")
drugs_filter <- c("AM-251", "AM-404", "AM-580", "aminoglutethimide", "aminopurvalanol-a")
LINCS_drug_signature_filenames <- overallDrugSignature$compute(gene_symbols, drugs_filter, drug_signature_cfg$perturbation_times, drug_signature_cfg$signature_base_dir, drug_signature_cfg$signature_base_filename)

LINCS_drug_signatures_ma <- overallMetanalysisDrugSignature6h24h$compute(LINCS_drug_signature_filenames)
LINCS_drug_signatures_ma_filename <- paste0(drug_signature_cfg$signature_base_dir, drug_signature_cfg$signature_ma_filename)
write.table(LINCS_drug_signatures_ma, paste0(LINCS_drug_signatures_ma_filename, ".csv"), sep = ";", row.names = F)
saveRDS(LINCS_drug_signatures_ma, paste0(LINCS_drug_signatures_ma_filename, ".Rds"))