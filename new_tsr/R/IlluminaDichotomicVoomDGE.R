IlluminaDichotomicVoomDGE <- R6Class(
  "IlluminaDichotomicVoomDGE",
  public = list(
    initialize = function(geneFilter = NA) {
      private$dichotomicVoomDGE <- DichotomicVoomDGE$new(geneFilter)
    },
    compute = function(rna_seq_filename, sample_01_map, test_sample_name, filter_by_protein_coding = F) {
      startTime <- Sys.time()
      dgrpLogger$log(sprintf("start %s differential gene expression computation", test_sample_name))
      rna_seq <- as.matrix(data.table::fread(rna_seq_filename, header = T, colClasses = "integer"), rownames = "GeneID")
      dge <- private$dichotomicVoomDGE$compute(rna_seq, sample_01_map, test_sample_name, filter_by_protein_coding)
      totalTime <- Sys.time() - startTime
      dgrpLogger$log(sprintf("end %s differential gene expression computation, time: %s %s", test_sample_name, totalTime, attr(totalTime, "units")))
      return(dge)
    }
  ),
  private = list(
    dichotomicVoomDGE = NA
  )
)
