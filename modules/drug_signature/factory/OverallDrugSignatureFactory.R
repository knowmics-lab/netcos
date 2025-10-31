library(R6)

source("modules/config.R")
source("modules/drug_signature/loader/GeneInfoListLoader.R")
source("modules/drug_signature/setuper/LINCSExperimentMetaDataSetuper.R")
source("modules/drug_signature/loader/RDS_LINCSExperimentDataLoader.R")
source("modules/drug_signature/julia/factory/JuliaDrugSignatureFactory.R")
source("modules/drug_signature/DrugSignature.R")
source("modules/drug_signature/DrugSignatureByScenarioSync.R")
source("modules/drug_signature/DrugSignatureByPerturbationTime.R")
source("modules/drug_signature/OverallDrugSignature.R")

OverallDrugSignatureFactory <- R6Class(
  "OverallDrugSignatureFactory",
  public = list(
    create = function(BLAS_num_threads = NA) {
      geneInfoListLoader <- GeneInfoListLoader$new(drug_signature_cfg$gene_info_filename)
      lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)
      rds_LINCSExperimentDataLoader <- RDS_LINCSExperimentDataLoader$new(drug_signature_cfg$lincs_splitted_level3_dir)
      juliaDrugSignatureFactory <- JuliaDrugSignatureFactory$new(rds_LINCSExperimentDataLoader)
      drugSignature <- juliaDrugSignatureFactory$create(BLAS_num_threads)
      drugSignatureByScenarioSync <- DrugSignatureByScenarioSync$new(drugSignature)
      drugSignatureByPerturbationTime <- DrugSignatureByPerturbationTime$new(lincsExperimentMetaDataSetuper, geneInfoListLoader, drugSignatureByScenarioSync)
      overallDrugSignature <- OverallDrugSignature$new(drugSignatureByPerturbationTime)
      return(overallDrugSignature)
    }
  )
)