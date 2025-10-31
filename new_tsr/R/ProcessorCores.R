
ProcessorCores <- R6Class(
  "ProcessorCores",
  public = list(
    initialize = function() {
      private$cores <- detectCores()
    },
    set = function(cores) {
      private$cores <- cores
    },
    get = function() {
      return(private$cores)
    },
    initCores = function(cores = NA) {
      if (!obj_is_na(cores)) {
        private$cores <- cores
      }
      cores <- processorCores$get()
      dgrpLogger$log(sprintf("processor cores: %s", cores))
      registerDoParallel(cores)
    }
  ),
  private = list(
    cores = NA
  )
)

processorCores <- ProcessorCores$new()