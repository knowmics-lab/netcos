library(R6)
library(readr)

source("modules/config.R")
source("modules/obj_is_na.R")

GeneInfoListLoader <- R6Class(
  "GeneInfoListLoader",
  public = list(
    initialize = function(gene_info_list_filename = NA) {
      if (obj_is_na(gene_info_list_filename))
        private$gene_info_list_filename <- drug_signature_cfg$gene_info_filename
      else
        private$gene_info_list_filename <- gene_info_list_filename
    },
    load = function() {
      gene_info_list <-
        read_delim(
          private$gene_info_list_filename,
          "\t",
          escape_double = FALSE,
          trim_ws = TRUE
        )
      return(gene_info_list)
    }
  ),
  private = list(
    gene_info_list_filename = NA
  )
)