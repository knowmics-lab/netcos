library(R6)
library(readr)
source("modules/TSRLogger.R")
RNASeqDataLoader <- R6Class(
  "RNASeqDataLoader",
  public = list(
    load = function(rna_seq_data_filename, rna_seq_metadata) {
      rna_seq_data <- as.matrix(data.table::fread(rna_seq_data_filename, header = T, colClasses = "integer"), rownames = "GeneID")
      rna_seq_data <- rna_seq_data[, rna_seq_metadata$sample]
      sizes <- dim(rna_seq_data)
      private$tsrLogger$log(sprintf("RNA Seq data loaded, size: %s X %s", sizes[1], sizes[2]))
      return(rna_seq_data)
    }
  ),
  private = list(
    tsrLogger = TSRLogger$new()
  )
)