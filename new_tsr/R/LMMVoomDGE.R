LMMVoomDGE <- R6Class(
  "LMMVoomDGE",
  public = list(
    initialize = function() {
      private$geneFilterByProteinCoding <- GeneFilterByProteinCoding$new()
      private$voomDGEMapper <- VoomDGEMapper$new()
    },
    compute = function(rna_seq_data, rna_seq_metadata, formula, filter_by_protein_coding = F) {
      data_size <- dim(rna_seq_data)
      dgrpLogger$log(sprintf("start LMM Voom differential gene expression computation, data size: %s X %s", data_size[1], data_size[2]))
      if (filter_by_protein_coding) {
        rna_seq_data <- private$geneFilterByProteinCoding$filterById(rna_seq_data)
        data_size <- dim(rna_seq_data)
        dgrpLogger$log(sprintf("RNA seq data size after filtering by protein coding genes: %s X %s", data_size[1], data_size[2]))
      }
      startTime <- Sys.time()
      dge <- DGEList(rna_seq_data, remove.zeros = TRUE)
      dge <- calcNormFactors(dge, method = 'upperquartile')
      parallelComputationParam <- SnowParam(processorCores$get(), "FORK", progressbar = TRUE)
      dgrpLogger$log("start voomWithDreamWeights computation")
      partialStartTime <- Sys.time()
      voom_data <- voomWithDreamWeights(dge, formula, rna_seq_metadata, BPPARAM = parallelComputationParam)
      totalTime <- Sys.time() - partialStartTime
      dgrpLogger$log(sprintf("end voomWithDreamWeights computation, time: %s %s", totalTime, attr(totalTime, "units")))
      dgrpLogger$log("start dream computation")
      partialStartTime <- Sys.time()
      differential_expression <- dream(voom_data, formula, rna_seq_metadata, BPPARAM = parallelComputationParam)
      totalTime <- Sys.time() - partialStartTime
      dgrpLogger$log(sprintf("end dream computation, time: %s %s", totalTime, attr(totalTime, "units")))
      differential_expression <- variancePartition::eBayes(differential_expression)
      differential_expression <- variancePartition::topTable(differential_expression, coef = 2, number = 10^6)
      differential_expression$std.error <- differential_expression$logFC / differential_expression$t
      totalTime <- Sys.time() - startTime
      dgrpLogger$log(sprintf("end LMM Voom differential gene expression computation, time: %s %s", totalTime, attr(totalTime, "units")))
      return(private$voomDGEMapper$map(differential_expression))
    }
  ),
  private = list(
    geneFilterByProteinCoding = NA,
    voomDGEMapper = NA
  )
)
