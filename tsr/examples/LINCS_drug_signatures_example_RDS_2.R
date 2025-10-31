source("modules/drug_signature/GeneSymbolsLoader.R")
source("modules/drug_signature/GeneInfoListLoader.R")
source("modules/drug_signature/LINCSExperimentMetaDataSetuper.R")
source("modules/drug_signature/RDS_LINCSExperimentDataLoader.R")
source("modules/drug_signature/DrugSignatureByScenarioAsync.R")
source("modules/drug_signature/DrugSignatureByPerturbationTime.R")
source("modules/drug_signature/OverallDrugSignature.R")

# setup
geneInfoListLoader <- GeneInfoListLoader$new(drug_signature_cfg$gene_info_filename)
lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)
rds_LINCSExperimentDataLoader <- RDS_LINCSExperimentDataLoader$new(drug_signature_cfg$lincs_splitted_level3_dir)
drugSignatureByScenarioAsync <- DrugSignatureByScenarioAsync$new(rds_LINCSExperimentDataLoader)
drugSignatureByPerturbationTime <- DrugSignatureByPerturbationTime$new(lincsExperimentMetaDataSetuper, geneInfoListLoader, drugSignatureByScenarioAsync)
overallDrugSignature <- OverallDrugSignature$new(drugSignatureByPerturbationTime)

print("start processing LINCS data")

gene_symbols <- c("DDR1", "PAX8", "FAU")
drugs_filter <- c("AM-251", "AM-404", "AM-580", "aminoglutethimide", "aminopurvalanol-a")
overallDrugSignature$compute(gene_symbols, drugs_filter, drug_signature_cfg$perturbation_times)