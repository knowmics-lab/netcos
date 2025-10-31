library(R6)
library(readr)

LINCSExperimentMetaDataLoader <- R6Class(
  "LINCSExperimentMetaDataLoader",
  public = list(
    initialize = function(experiments_meta_data_filename) {
      private$experiments_meta_data_filename <- experiments_meta_data_filename
    },
    load = function() {
      experimente_meta_data <-
        read_delim(
          private$experiments_meta_data_filename,
          "\t",
          escape_double = FALSE,
          trim_ws = TRUE
        )
      return(experimente_meta_data[, c("inst_id", "rna_plate", "pert_id", "pert_iname", "pert_type", "pert_dose", "pert_dose_unit", "pert_time", "pert_time_unit", "cell_id")])
    }
  ),
  private = list(
    experiments_meta_data_filename = NA
  )
)