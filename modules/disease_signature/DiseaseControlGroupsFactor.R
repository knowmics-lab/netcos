library(R6)

DiseaseControlGroupsFactor <- R6Class(
    "DiseaseControlGroupsFactor",
    public = list(
      create = function(samples_groups_map_filter, disease_name) {
        disease_control_map <- subset(samples_groups_map_filter, samples_groups_map_filter != "X")
        disease_control_map[which(disease_control_map == '0')] <- disease_name
        disease_control_map[which(disease_control_map == '1')] <- "Control"
        disease_control_map <- factor(disease_control_map, levels = c("Control", disease_name))
        return(disease_control_map)
      }
    )
  )