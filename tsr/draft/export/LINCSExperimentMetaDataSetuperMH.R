library(R6)
source("modules/obj_is_na.R")
source("modules/config.R")
source("modules/drug_signature/LINCSExperimentMetaDataLoader.R")

LINCSExperimentMetaDataSetuperMH <- R6Class(
  "LINCSExperimentMetaDataSetuperMH",
  public = list(
    setup = function(perturbation_times, drugs_filter = NA) {
      # "data/LINCS/GSE92742_Broad_LINCS_inst_info.txt"
      # experiment = instance
      raw_experiments_meta_data <- private$lincsExperimentMetaDataLoader$load()

      # la condizione pert_dose_unit != "-666" è superflua perchè non ci sono esperimenti
      # per cui vale pert_dose='10' and pert_dose_unit='-666'
      experiments_meta_data <- subset(raw_experiments_meta_data,
                                      pert_type == "trt_cp" &
                                        pert_dose == 10 &
                                        pert_time %in% perturbation_times
      )

      # calcola la frequenza di linee cellulari per farmaco
      cell_lines_per_drug_freq <- as.data.frame(table(unique(experiments_meta_data[, c("cell_id", "pert_iname")])$pert_iname), stringsAsFactors = FALSE)
      # elimina le linee cellulari che appaiono con una frequenza minore di 5
      cell_lines_per_drug_freq <- cell_lines_per_drug_freq[cell_lines_per_drug_freq$Freq >= 5,]

      drugs <- unique(cell_lines_per_drug_freq$Var1)

      if (!obj_is_na(drugs_filter)) {
        drugs <- drugs[drugs %in% drugs_filter]
      }
      # elimina gli esperimenti sui quali non sono testati i farmaci in drug_names
      experiments_meta_data <- experiments_meta_data[experiments_meta_data$pert_iname %in% drugs,]
      # estrae tutti gli esperimenti eseguiti sulle stesse piastre in filtered_experiment_info_list
      raw_experiments_meta_data <- raw_experiments_meta_data[raw_experiments_meta_data$rna_plate %in% unique(experiments_meta_data$rna_plate),]

      # aggiunge il veicolo di controllo
      experiments_meta_data <- rbind(experiments_meta_data, raw_experiments_meta_data[raw_experiments_meta_data$pert_type == "ctl_vehicle" & raw_experiments_meta_data$pert_time %in% perturbation_times,])

      # le prossime due righe ci assicurano che il DMSO sia in testa al factor
      # l'elemento in testa al factor non appare nell'elenco "drug_coefs"
      # il primo elemento del factor viene considerato come livello di espressione genica
      # di riferimento per il calcolo del log2 fold change
      drugs <- c("DMSO", drugs[drugs != "DMSO"])
      experiments_meta_data$pert_iname <- factor(experiments_meta_data$pert_iname, levels = drugs)
      #drug_names <- unique(experiments_meta_data$pert_iname)
      #drug_names <- c("DMSO",drug_names[drug_names != "DMSO"])
      #experiments_meta_data$pert_iname <- factor(experiments_meta_data$pert_iname, levels = drug_names)
      experiments_meta_data <- experiments_meta_data[order(experiments_meta_data$inst_id),]
      experiments_env <- new.env()
      experiments_env$experiments_meta_data <- experiments_meta_data[, c("inst_id", "rna_plate", "rna_well", "pert_id", "pert_iname", "pert_type", "pert_dose", "pert_dose_unit", "pert_time", "pert_time_unit", "cell_id")]
      experiments_env$experiments_meta_data <- experiments_env$experiments_meta_data[order(experiments_env$experiments_meta_data$inst_id),]
      return(experiments_env)
    }
  ),
  private = list(
    lincsExperimentMetaDataLoader = LINCSExperimentMetaDataLoader$new(drug_signature_cfg$experiments_meta_data_filename)
  )
)