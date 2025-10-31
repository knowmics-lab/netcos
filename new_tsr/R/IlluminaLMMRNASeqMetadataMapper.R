
IlluminaLMMRNASeqMetadataMapper <- R6Class(
  "IlluminaLMMRNASeqMetadataMapper",
  public = list(
    map = function(rna_seq_metadata, tissue_statuses_to_be_tested, tissue_statuses_map) {
      rna_seq_metadata <- rna_seq_metadata[, c("accession", "tissue", "tissue_status")]
      rna_seq_metadata <- subset(rna_seq_metadata, tissue_status %in% tissue_statuses_to_be_tested)
      for (i in 1:length(tissue_statuses_to_be_tested)) {
        rna_seq_metadata[rna_seq_metadata == tissue_statuses_to_be_tested[i]] <- tissue_statuses_map[i]
      }
      colnames(rna_seq_metadata)[1] <- "sample"
      rna_seq_metadata$tissue_status <- factor(rna_seq_metadata$tissue_status, levels = tissue_statuses_map)
      return(rna_seq_metadata)
    }
  )
)