library(R6)
source("modules/config.R")
source("modules/obj_is_na.R")

DrugSignatureInputArguments <- R6Class(
  "DrugSignatureInputArguments",
  public = list(
    initialize = function(default_genes_filename = NA) {
      private$default_genes_filename <- default_genes_filename
    },
    get = function() {
      input_args <- commandArgs(trailingOnly = TRUE)
      if (!obj_is_na(input_args[1])) {
        genes_filename <- input_args[1]
      }else {
        genes_filename <- private$default_genes_filename
      }
      if (!obj_is_na(input_args[2])) {
        perturbation_times <- unlist(strsplit(input_args[2], ","))
      }else {
        perturbation_times <- drug_signature_cfg$perturbation_times
      }
      if (!obj_is_na(input_args[3])) {
        BLAS_num_threads <- unlist(strsplit(input_args[3], ","))
      }else {
        BLAS_num_threads <- tsr_julia_cfg$BLAS_num_threads
      }
      if (!obj_is_na(input_args[4])) {
        num_chunks <- strtoi(input_args[4])
      }else {
        num_chunks <- parallel_computation$max_cores
      }
      return(list(
        genes_filename = genes_filename,
        perturbation_times = perturbation_times,
        BLAS_num_threads = BLAS_num_threads,
        num_chunks = num_chunks
      ))
    }
  ),
  private = list(
    default_genes_filename = NA
  )
)