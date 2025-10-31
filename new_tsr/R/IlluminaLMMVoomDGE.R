IlluminaLMMVoomDGE <- R6Class(
  "IlluminaLMMVoomDGE",
  public = list(
    initialize = function(geneFilter = NA) {
      private$illuminaLMMRNASeqLoader <- IlluminaLMMRNASeqLoader$new(geneFilter)
      private$lmmVoomDGE <- LMMVoomDGE$new()
    },
    compute = function(rna_seq_metadata_filename, rna_seq_data_filename, tissue_statuses_to_be_tested, tissue_statuses_map, formula, filter_by_protein_coding = F) {
      startTime <- Sys.time()
      dgrpLogger$log(sprintf("start differential gene expression computation"))
      rna_seq <- private$illuminaLMMRNASeqLoader$load(rna_seq_metadata_filename, rna_seq_data_filename, tissue_statuses_to_be_tested, tissue_statuses_map)
      dge <- private$lmmVoomDGE$compute(rna_seq$data, rna_seq$metadata, formula, filter_by_protein_coding)
      totalTime <- Sys.time() - startTime
      dgrpLogger$log(sprintf("end differential gene expression computation, time: %s %s", totalTime, attr(totalTime, "units")))
      return(dge)
    }
  ),
  private = list(
    illuminaLMMRNASeqLoader = NA,
    lmmVoomDGE = NA
  )
)