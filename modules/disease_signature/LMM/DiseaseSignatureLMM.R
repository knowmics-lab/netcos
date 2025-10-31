library(R6)
library(edgeR)
library(variancePartition)

source("modules/config.R")
source("modules/TSRLogger.R")
source("modules/disease_signature/GeneFilterByProteinCoding.R")

DiseaseSignatureLMM <- R6Class(
  "DiseaseSignatureLMM",
  public = list(
    compute = function(rna_seq_data, rna_seq_metadata, formula, filter_by_proteing_coding = F) {
      if (filter_by_proteing_coding) {
        data_size <- dim(rna_seq_data)
        private$tsrLogger$log(sprintf("LMM disease signature data size before filter protein coding genes: %s X %s", data_size[1], data_size[2]))
        rna_seq_data <- private$geneFilterByProteinCoding$filter(rna_seq_data)
      }
      data_size <- dim(rna_seq_data)
      private$tsrLogger$log(sprintf("start LMM disease signature computation, data size: %s X %s", data_size[1], data_size[2]))
      startTime <- Sys.time()
      dge <- DGEList(rna_seq_data, remove.zeros = TRUE)
      dge <- calcNormFactors(dge, method = 'upperquartile')
      parallelComputationParam <- SnowParam(parallel_computation$LMM_disease_signature_max_cores, "FORK", progressbar = TRUE)
      private$tsrLogger$log("start voomWithDreamWeights computation")
      partialStartTime <- Sys.time()
      voom_data <- voomWithDreamWeights(dge, formula, rna_seq_metadata, BPPARAM = parallelComputationParam)
      totalTime <- Sys.time() - partialStartTime
      private$tsrLogger$log(sprintf("end voomWithDreamWeights computation, time: %s %s", totalTime, attr(totalTime, "units")))
      private$tsrLogger$log("start dream computation")
      partialStartTime <- Sys.time()
      differential_expression <- dream(voom_data, formula, rna_seq_metadata, BPPARAM = parallelComputationParam)
      totalTime <- Sys.time() - partialStartTime
      private$tsrLogger$log(sprintf("end dream computation, time: %s %s", totalTime, attr(totalTime, "units")))
      differential_expression <- variancePartition::eBayes(differential_expression)
      disease_signature <- variancePartition::topTable(differential_expression, coef = 2, number = 10^6)
      totalTime <- Sys.time() - startTime
      private$tsrLogger$log(sprintf("end LMM disease signature computation, time: %s %s", totalTime, attr(totalTime, "units")))
      return(disease_signature)
    }
  ),
  private = list(
    tsrLogger = TSRLogger$new(),
    geneFilterByProteinCoding = GeneFilterByProteinCoding$new()
  )
)