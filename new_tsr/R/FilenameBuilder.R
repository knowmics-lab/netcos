
FilenameBuilder <- R6Class(
  "FilenameBuilder",
  public = list(
    valid_filename_pattern = function(filename_pattern) {
      return(grepl("#id", filename_pattern))
    },
    build = function(filename_pattern, gene_id, gene_symbol = NA) {
      filename <- filename_pattern
      if (!is.na(gene_symbol))
        filename <- gsub("#symbol", gene_symbol, filename)
      if (!is.na(gene_id))
        filename <- gsub("#id", gene_id, filename)
      return(filename)
    }
  )
)