library(R6)
DrugSignatureInputChunksBuilder <- R6Class(
  "DrugSignatureInputChunksBuilder",
  public = list(
    build = function(gene_sinbols_list, perturbation_times, BLAS_num_threads_list) {
      chunks <- list()
      gene_sinbols_list_size <- length(gene_sinbols_list)
      perturbation_times_size <- length(perturbation_times)
      for (i in 1:perturbation_times_size) {
        for (j in 1:gene_sinbols_list_size) {
          pos <- (i - 1) * gene_sinbols_list_size + j
          chunk <- list(
            gene_symbols = unlist(gene_sinbols_list[j]),
            perturbation_time = perturbation_times[i],
            BLAS_num_threads = BLAS_num_threads_list[pos],
            number = pos
          )
          chunks[[pos]] <- chunk
        }
      }
      return(chunks)
    }
  )
)