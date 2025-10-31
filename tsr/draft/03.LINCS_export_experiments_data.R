source("modules/config.R")
source("modules/drug_signature/loader/GeneSymbolsLoader.R")
source("modules/drug_signature/loader/GeneInfoListLoader.R")
source("modules/drug_signature/setuper/LINCSExperimentMetaDataSetuper.R")
source("modules/drug_signature/loader/RDS_LINCSExperimentDataLoader.R")
source("modules/drug_signature/DrugSignatureExportExperimentData.R")
source("modules/drug_signature/DrugSignatureByScenarioAsync.R")
source("modules/drug_signature/DrugSignatureByPerturbationTime.R")
source("modules/drug_signature/OverallDrugSignature.R")
source("modules/drug_signature/DrugSignatureInputArguments.R")
source("modules/TSRLogger.R")

# setup
geneSymbolsLoader <- GeneSymbolsLoader$new()
geneInfoListLoader <- GeneInfoListLoader$new(drug_signature_cfg$gene_info_filename)
lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)
rds_LINCSExperimentDataLoader <- RDS_LINCSExperimentDataLoader$new(drug_signature_cfg$lincs_splitted_level3_dir)

#drug_signature_cfg$csv_experiments_data <- "output/drug_signature/LINCS/experiments_data/"
drugSignature <- DrugSignatureExportExperimentData$new(rds_LINCSExperimentDataLoader, drug_signature_cfg$csv_experiments_data)
drugSignatureByScenarioSync <- DrugSignatureByScenarioAsync$new(drugSignature)
drugSignatureByPerturbationTime <- DrugSignatureByPerturbationTime$new(lincsExperimentMetaDataSetuper, geneInfoListLoader, drugSignatureByScenarioSync)
overallDrugSignature <- OverallDrugSignature$new(drugSignatureByPerturbationTime)
drugSignatureInputArguments <- DrugSignatureInputArguments$new()
tsrLogger <- TSRLogger$new()

# computation
tsrLogger$log("start processing LINCS data")
inputArguments <- drugSignatureInputArguments$get()
tsrLogger$log(sprintf("gene list file: %s", inputArguments$genes_filename))
perturbation_times_label <- paste0(paste(inputArguments$perturbation_times, collapse = "h, "), "h")
tsrLogger$log(sprintf("perturbation times: %s", perturbation_times_label))
gene_symbols <- geneSymbolsLoader$load(inputArguments$genes_filename)
drugs_filter <- read.table("data/drugs.csv", header = T)$drug
overallDrugSignature$compute(gene_symbols, drugs_filter, inputArguments$perturbation_times)