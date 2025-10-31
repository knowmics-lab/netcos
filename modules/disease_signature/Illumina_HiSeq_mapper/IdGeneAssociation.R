library(R6)
source("modules/config.R")

IdGeneAssociation <- R6Class(
  "IdGeneAssociation",
  public = list(
    load = function() {
      id_gene_association <- data.table::fread(private$annot_filename, header = T, quote = "", stringsAsFactors = F, data.table = F)
      colnames(id_gene_association)[1] <- "ID"
      colnames(id_gene_association)[2] <- "gene"
      id_gene_association <- id_gene_association[, c("ID", "gene")]
      return(id_gene_association)
    }
  ),
  private = list(
    annot_filename = disease_signature_cfg$Illumina_HiSeq_annot_filename
  )
)