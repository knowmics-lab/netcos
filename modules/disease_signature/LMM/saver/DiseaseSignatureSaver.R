library(R6)

DiseaseSignatureSaver <- R6Class(
  "DiseaseSignatureSaver",
  public = list(
    save = function(disease_signature, filename) {
      write.table(disease_signature, file = filename, sep = ";", row.names = FALSE, dec = ",")
      saveRDS(disease_signature, gsub(".csv", ".Rds", filename))
    }
  )
)