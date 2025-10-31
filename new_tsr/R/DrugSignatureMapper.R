
DrugSignatureMapper <- R6Class(
  "DrugSignatureMapper",
  public = list(
    map = function(drug_signature, t_value_column_name) {
      rownames(drug_signature) <- drug_signature$gene_id
      drug_signature <- drug_signature[, t_value_column_name, drop = F]
      colnames(drug_signature) <- "estimate"
      return(drug_signature)
    }
  )
)