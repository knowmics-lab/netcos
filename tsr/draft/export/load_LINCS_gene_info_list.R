library(readr)

load_LINCS_gene_info_list <- function(gene_info_list_filename) {
  gene_info_list <-
    read_delim(
      gene_info_list_filename,
      "\t",
      escape_double = FALSE,
      trim_ws = TRUE
    )
  return(gene_info_list)
}