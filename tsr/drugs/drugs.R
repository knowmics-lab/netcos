source("modules/config.R")
source("modules/drug_signature/LINCSExperimentMetaDataSetuper.R")

lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)
experiments_6 <- lincsExperimentMetaDataSetuper$setup("6")
experiments_24 <- lincsExperimentMetaDataSetuper$setup("24")
drugs_6 <- unique(experiments_6$experiments_meta_data$pert_iname)
drugs_24 <- unique(experiments_24$experiments_meta_data$pert_iname)
drugs <- intersect(drugs_6, drugs_24)
drugs <- subset(drugs, drugs != "DMSO")
drugs <- drugs[order(drugs)]

write.table(drugs, "data/drugs.csv", quote = F, row.names = F, col.names = "drug")
print("end computation")