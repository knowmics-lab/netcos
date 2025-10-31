library(R6)
source("modules/disease_signature/DiseaseControlGroupsFactor.R")

SamplesGroups <- R6Class(
  "SamplesGroups",
  public = list(
    groups = function(samples_groups_map, disease_name) {
      samples_groups_map_filter <- strsplit(samples_groups_map, split = "")[[1]]
      disease_control_groups <- private$diseaseControlGroupsFactor$create(samples_groups_map_filter, disease_name)
      sample_groups_idxes <- which(samples_groups_map_filter != "X")
      return(
        list(
          disease_control_groups = disease_control_groups,
          sample_groups_idxes = sample_groups_idxes
        )
      )
    }
  ),
  private = list(
    diseaseControlGroupsFactor = DiseaseControlGroupsFactor$new()
  )
)

