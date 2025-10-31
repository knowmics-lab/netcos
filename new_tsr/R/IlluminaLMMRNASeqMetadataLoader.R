
IlluminaLMMRNASeqMetadataLoader <- R6Class(
  "IlluminaLMMRNASeqMetadataLoader",
  public = list(
    initialize = function() {
      private$IlluminaLMMRNASeqMetadataMapper <- IlluminaLMMRNASeqMetadataMapper$new()
    },
    load = function(rna_seq_metadata_filename, tissue_statuses_to_be_tested, tissue_statuses_map) {
      rna_seq_metadata <- read_delim(rna_seq_metadata_filename, delim = ";")
      sizes <- dim(rna_seq_metadata)
      dgrpLogger$log(sprintf("RNA Seq metadata loaded, size: %s X %s", sizes[1], sizes[2]))
      rna_seq_metadata <- private$IlluminaLMMRNASeqMetadataMapper$map(rna_seq_metadata, tissue_statuses_to_be_tested, tissue_statuses_map)
      return(rna_seq_metadata)
    }
  ),
  private = list(
    IlluminaLMMRNASeqMetadataMapper = NA
  )
)