
DrugSignatureLoaderByDrugName <- R6Class(
  "DrugSignatureLoaderByDrugName",
  inherit = DrugSignatureLoaderAbstract,
  public = list(
    initialize = function(drug_dge_dir, t_value_column_name = "t.value", disease_drug_common_genes = NA) {
      private$drug_dge_dir <- drug_dge_dir
      private$t_value_column_name <- t_value_column_name
      private$disease_drug_common_genes <- disease_drug_common_genes
      private$drugSignatureMapper <- DrugSignatureMapper$new()
    },
    load = function(drug_name) {
      filename <- paste0(private$drug_dge_dir, drug_name, ".Rds")
      drug_dge <- readRDS(filename)
      if (!obj_is_na(private$disease_drug_common_genes)) {
        drug_signature <- subset(drug_dge, drug_dge$gene_id %in% private$disease_drug_common_genes)
      }else {
        drug_signature <- drug_dge
      }
      return(private$drugSignatureMapper$map(drug_signature, private$t_value_column_name))
    }
  ),
  private = list(
    drug_dge_dir = NA,
    t_value_column_name = NA,
    drugSignatureMapper = NA
  )
)