
LINCSMetadataSetuper <- R6Class(
  "LINCSMetadataSetuper",
  inherit = LINCSMetadataSetuperAbstract,
  public = list(
    setup = function(perturbation_times, drugs_filter = NA) {
      # "data/LINCS/GSE92742_Broad_LINCS_inst_info.txt"
      # experiment = instance
      raw_metadata <- package_readRDS(config$LINCS_metadata_RDS_filename)[, c("inst_id", "rna_plate", "pert_id", "pert_iname", "pert_type", "pert_dose", "pert_time", "cell_id")]

      # la condizione pert_dose_unit != "-666" è superflua perchè non ci sono esperimenti
      # per cui vale pert_dose='10' and pert_dose_unit='-666'
      metadata <- subset(raw_metadata,
                         pert_type == "trt_cp" &
                           pert_dose == 10 &
                           pert_time %in% perturbation_times
      )

      # calcola la frequenza di linee cellulari per farmaco
      cell_lines_per_drug_freq <- as.data.frame(table(unique(metadata[, c("cell_id", "pert_iname")])$pert_iname), stringsAsFactors = FALSE)
      # elimina le linee cellulari che appaiono con una frequenza minore di 5
      cell_lines_per_drug_freq <- cell_lines_per_drug_freq[cell_lines_per_drug_freq$Freq >= 5,]

      drugs <- unique(cell_lines_per_drug_freq$Var1)

      if (!obj_is_na(drugs_filter)) {
        drugs <- drugs[drugs %in% drugs_filter]
      }
      # elimina gli esperimenti sui quali non sono testati i farmaci in drug_names
      metadata <- metadata[metadata$pert_iname %in% drugs,]
      # estrae tutti gli esperimenti eseguiti sulle stesse piastre in metadata
      raw_metadata <- raw_metadata[raw_metadata$rna_plate %in% unique(metadata$rna_plate),]

      # aggiunge il veicolo di controllo
      metadata <- rbind(metadata, raw_metadata[raw_metadata$pert_type == "ctl_vehicle" & raw_metadata$pert_time %in% perturbation_times,])

      # le prossime due righe ci assicurano che il DMSO sia in testa al factor
      # l'elemento in testa al factor non appare nell'elenco "drug_coefs"
      # il primo elemento del factor viene considerato come livello di espressione genica
      # di riferimento per il calcolo del log2 fold change
      drugs <- c("DMSO", drugs[drugs != "DMSO"])
      metadata$pert_iname <- factor(metadata$pert_iname, levels = drugs)
      #drug_names <- unique(metadata$pert_iname)
      #drug_names <- c("DMSO",drug_names[drug_names != "DMSO"])
      #metadata$pert_iname <- factor(metadata$pert_iname, levels = drug_names)
      metadata <- metadata[, c("inst_id", "rna_plate", "pert_id", "pert_iname", "pert_time", "cell_id")]
      metadata <- metadata[order(metadata$inst_id),]
      return(metadata)
    }
  )
)