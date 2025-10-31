
RandomConnectivityScoreDistributionParallel <- R6Class(
  "RandomConnectivityScoreDistributionParallel",
  inherit = RandomConnectivityScoreDistributionAbstract,
  public = list(
    initialize = function() {
      private$randomConnectivityScoreDistribution <- RandomConnectivityScoreDistribution$new()
    },
    compute = function(n_disease_up_regulated_genes, n_disease_down_regulated_genes, n_drug_signatures_genes, n_permutations) {
      cores <- processorCores$get()
      startTime <- Sys.time()
      dgrpLogger$log("start random connectivity score computation")
      block_size <- floor(n_permutations / cores)
      first_block_size <- n_permutations - block_size * (cores - 1)
      block_sizes <- rep(block_size, cores)
      block_sizes[1] <- first_block_size

      output <- foreach(i = 1:cores, .combine = c) %dopar% {
        private$
          randomConnectivityScoreDistribution$
          compute(n_disease_up_regulated_genes, n_disease_down_regulated_genes, n_drug_signatures_genes, block_sizes[i])
      }
      totalTime <- Sys.time() - startTime
      dgrpLogger$log(sprintf("end random connectivity score computation, time: %s %s", totalTime, attr(totalTime, "units")))
      return(output)
    }
  ),
  private = list(
    randomConnectivityScoreDistribution = NA
  )
)
