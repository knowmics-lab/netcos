source("modules/config.R")
source("modules/TSRLogger.R")
source("modules/disease_signature/Illumina_HiSeq_mapper/IlluminaHiSeqMapper.R")
source("modules/disease_signature/DiseaseSignature.R")
source("modules/disease_signature/filter/LowCountsGeneFilter.R")

# setup
# download genetic db from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE3nnn/GSE3307/matrix/
tsrLogger <- TSRLogger$new()
lowCountsGeneFilter <- LowCountsGeneFilter$new()
illuminaHiSeqMapper <- IlluminaHiSeqMapper$new(lowCountsGeneFilter)
diseaseSignature <- DiseaseSignature$new()

# computation
startTime <- Sys.time()
tsrLogger$log("start ipf disease signature computation")

ipf_gene_experiments_data <- illuminaHiSeqMapper$map(ipf_disease_signature_cfg$gse_92592_filename, ipf_disease_signature_cfg$samples_groups_map, "ipf")
ipf_signature <- diseaseSignature$compute(ipf_gene_experiments_data, F)
filename_no_ext <- paste0(ipf_disease_signature_cfg$signature_base_dir, ipf_disease_signature_cfg$signature_filename)
write.table(ipf_signature, file = paste0(filename_no_ext, ".csv"), sep = ";", row.names = FALSE, dec = ",")
saveRDS(ipf_signature, paste0(filename_no_ext, ".Rds"))

totalTime <- Sys.time() - startTime
tsrLogger$log(sprintf("end ipf disease signature computation, time: %s %s", totalTime, attr(totalTime, "units")))

