library(R6)
library(foreach)

source("modules/config.R")
source("modules/TSRLogger.R")
source("modules/CoreProcessorInit.R")
source("modules/connectivity_score/cmap_score_new/RandomConnectivityScore.R")

RandomConnectivityScoreParallel <- R6Class(
  "RandomConnectivityScoreParallel",
  public = list(
    compute = function(n_genes_up, n_genes_down, n_genes, n_permutations) {
      startTime <- Sys.time()
      private$tsrLogger$log("start random connectivity score computation")
      cores <- parallel_computation$max_cores
      block_size <- floor(n_permutations / cores)
      first_block_size <- n_permutations - block_size * (cores - 1)
      block_sizes <- rep(block_size, cores)
      block_sizes[1] <- first_block_size

      output <- foreach(i = 1:cores, .combine = c) %dopar% {
        private$randomConnectivityScore$compute(n_genes_up, n_genes_down, n_genes, block_sizes[i])
      }
      totalTime <- Sys.time() - startTime
      private$tsrLogger$log(sprintf("end random connectivity score computation, time: %s %s", totalTime, attr(totalTime, "units")))
      return(output)
    }
  ), private = list(
    tsrLogger = TSRLogger$new(),
    randomConnectivityScore = RandomConnectivityScore$new(),
    coreProcessorInit = CoreProcessorInit$new()
  )
)
