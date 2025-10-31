library(R6)
library(doParallel)
source("modules/config.R")
source("modules/TSRLogger.R")

CoreProcessorInit <- R6Class(
  public = list(
    init = function() {
      parallel_computation$max_cores
      coreNumbers <- detectCores()
      private$tsrLogger$log(sprintf("processor core numbers: %s", coreNumbers))
      usedCores <- min(coreNumbers, parallel_computation$max_cores)
      private$tsrLogger$log(sprintf("processor cores used: %s", usedCores))
      registerDoParallel(usedCores)
    }
  ),
  private = list(
    tsrLogger = TSRLogger$new()
  )
)