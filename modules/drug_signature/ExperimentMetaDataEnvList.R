library(R6)
source("modules/obj_is_na.R")

ExperimentMetaDataEnvList <- R6Class(
  "ExperimentMetaDataEnvList",
  public = list(
    initialize = function(experimentMetaDataSetuper) {
      if (!"ExperimentMetaDataSetuper" %in% class(experimentMetaDataSetuper))
        stop("the experimentMetaDataSetuper class type must be a subclass of ExperimentsMetaDataSetuper")
      private$experimentMetaDataSetuper <- experimentMetaDataSetuper
    },
    get = function(perturbation_times, drugs_filter, cell_name) {
      experiments_env_list <- list()
      for (perturbation_time in perturbation_times) {
        experiments_env_list[[perturbation_time]] <- private$experimentMetaDataSetuper$setup(perturbation_time, drugs_filter,cell_name)
      }
      return(experiments_env_list)
    }
  ),
  private = list(
    experimentMetaDataSetuper = NA
  )
)
