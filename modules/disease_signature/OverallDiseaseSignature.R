library(R6)
source("modules/TSRLogger.R")
source("modules/obj_is_na.R")
source("modules/disease_signature/Illumina_HiSeq_mapper/IlluminaHiSeqMapper.R")
source("modules/disease_signature/filter/LowCountsGeneFilter.R")
source("modules/disease_signature/DiseaseSignature.R")

OverallDiseaseSignature <- R6Class(
  "OverallDiseaseSignature",
  public = list(
    initialize = function(geneFilter = NA) {
      private$tsrLogger <- TSRLogger$new()
      private$diseaseSignature <- DiseaseSignature$new()
      if (obj_is_na(geneFilter)) {
        geneFilter <- LowCountsGeneFilter$new()
      }
      private$illuminaHiSeqMapper <- IlluminaHiSeqMapper$new(geneFilter)
    },
    compute = function(rna_seq_filename, samples_groups_map, disease_name, signature_filename) {
      startTime <- Sys.time()
      private$tsrLogger$log(sprintf("start %s signature computation", disease_name))
      gene_experiments_data <- private$illuminaHiSeqMapper$map(rna_seq_filename, samples_groups_map, disease_name)
      signature <- private$diseaseSignature$compute(gene_experiments_data, F)
      write.table(signature, file = signature_filename, sep = ";", row.names = FALSE, dec = ",")
      saveRDS(signature, gsub(".csv", ".Rds", signature_filename))
      totalTime <- Sys.time() - startTime
      private$tsrLogger$log(sprintf("end %s signature computation, time: %s %s", disease_name, totalTime, attr(totalTime, "units")))
    }
  ),
  private = list(
    tsrLogger = NA,
    diseaseSignature = NA,
    illuminaHiSeqMapper = NA
  )
)