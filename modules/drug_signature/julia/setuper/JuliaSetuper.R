library(R6)
library(JuliaCall)
source("modules/config.R")
source("modules/obj_is_na.R")
source("modules/TSRLogger.R")

JuliaSetuper <- R6Class(
  "JuliaSetuper",
  public = list(
    setup = function(BLAS_num_threads = NA) {
      private$tsrLogger$log("julia setup...")
      julia_setup()
      julia_library("MixedModels")
      if (obj_is_na(BLAS_num_threads)) {
        BLAS_num_threads <- tsr_julia_cfg$BLAS_num_threads
      }
      julia_command(sprintf("BLAS.set_num_threads(%s)", BLAS_num_threads))
      private$tsrLogger$log(paste0("set OPENBLAS_NUM_THREADS = ", BLAS_num_threads))
      private$tsrLogger$log(paste0("check OPENBLAS_NUM_THREADS = ", julia_eval("BLAS.get_num_threads()")))

      return(NA)
    }
  ),
  private = list(
    tsrLogger = TSRLogger$new()
  )
)