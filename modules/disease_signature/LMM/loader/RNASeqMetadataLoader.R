library(R6)
library(readr)
source("modules/TSRLogger.R")

RNASeqMetadataLoader <- R6Class(
  "RNASeqMetadataLoader",
  public = list(
    load = function(rna_seq_metadata_filename) {
      rna_seq_metadata <- read_delim(rna_seq_metadata_filename, delim = ";")
      sizes <- dim(rna_seq_metadata)
      private$tsrLogger$log(sprintf("RNA Seq metadata loaded, size: %s X %s", sizes[1], sizes[2]))
      return(rna_seq_metadata)
    }
  ),
  private = list(
    tsrLogger = TSRLogger$new()
  )
)