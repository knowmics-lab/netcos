source("modules/config.R")
source("modules/TSRLogger.R")
source("modules/disease_signature/als_1/ALS_GSE3307_GPL97DataPrepare.R")
source("modules/disease_signature/DiseaseSignature.R")

# setup
# download genetic db from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE3nnn/GSE3307/matrix/
tsrLogger <- TSRLogger$new()
diseaseSignature <- DiseaseSignature$new()
als_GSE3307_GPL97DataPrepare <- ALS_GSE3307_GPL97DataPrepare$new()

# computation
startTime <- Sys.time()
tsrLogger$log("start als 1 disease signature computation")

als_1_gene_experiments_data <- als_GSE3307_GPL97DataPrepare$prepare(als_1_disease_signature_cfg$gse_3307_GPL97_filename, als_1_disease_signature_cfg$gse_3307_GPL97_annot_dir, als_1_disease_signature_cfg$samples_groups_map)
als_1_signature <- diseaseSignature$compute(als_1_gene_experiments_data, F)
filename_no_ext <- paste0(als_1_disease_signature_cfg$signature_base_dir, als_1_disease_signature_cfg$signature_filename)
write.table(als_1_signature, file = paste0(filename_no_ext, ".csv"), sep = ";", row.names = FALSE, dec = ",")
saveRDS(als_1_signature, paste0(filename_no_ext, ".Rds"))

totalTime <- Sys.time() - startTime
tsrLogger$log(sprintf("end als 1 disease signature computation, time: %s %s",totalTime, attr(totalTime, "units")))
