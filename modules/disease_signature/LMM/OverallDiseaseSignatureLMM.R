library(R6)
source("modules/TSRLogger.R")
source("modules/disease_signature/LMM/loader/RNASeqLoader.R")
source("modules/disease_signature/LMM/DiseaseSignatureLMM.R")
source("modules/disease_signature/LMM/mapper/DiseaseSignatureMapper.R")
source("modules/disease_signature/LMM/saver/DiseaseSignatureSaver.R")

OverallDiseaseSignatureLMM <- R6Class(
  "OverallDiseaseSignatureLMM",
  public = list(
    initialize = function(geneFilter) {
      private$tsrLogger <- TSRLogger$new()
      private$rnaSeqLoader <- RNASeqLoader$new(geneFilter)
      private$diseaseSignatureLMM <- DiseaseSignatureLMM$new()
      private$diseaseSignatureMapper <- DiseaseSignatureMapper$new()
      private$diseaseSignatureSaver <- DiseaseSignatureSaver$new()
    },
    compute = function(rna_seq_metadata_filename, rna_seq_data_filename, tissue_statuses_to_be_tested, tissue_statuses_map, formula, signature_filename) {
      startTime <- Sys.time()
      private$tsrLogger$log(sprintf("start disease signature computation"))
      rna_seq <- private$rnaSeqLoader$load(rna_seq_metadata_filename, rna_seq_data_filename, tissue_statuses_to_be_tested, tissue_statuses_map)
      disease_signature <- private$diseaseSignatureLMM$compute(rna_seq$data, rna_seq$metadata, formula)
      disease_signature <- private$diseaseSignatureMapper$map(disease_signature)
      private$diseaseSignatureSaver$save(disease_signature, signature_filename)
      totalTime <- Sys.time() - startTime
      private$tsrLogger$log(sprintf("end disease signature computation, time: %s %s", totalTime, attr(totalTime, "units")))
    }
  ),
  private = list(
    tsrLogger = NA,
    rnaSeqLoader = NA,
    diseaseSignatureLMM = NA,
    diseaseSignatureMapper = NA,
    diseaseSignatureSaver = NA
  )
)