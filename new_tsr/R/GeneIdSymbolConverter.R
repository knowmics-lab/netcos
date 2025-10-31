
GeneIdSymbolConverter <- R6Class(
  "GeneIdSymbolConverter",
  public = list(
    initialize = function(gene_set = config$gene_set_Human_GRCh38.p13) {
      private$gene_set <- package_readRDS(paste0(config$gene_set_base_path, gene_set, config$gene_set_filename_suffix))
    },

    idToSymbol = function(gene_id_list) {
      return(private$fromFiled1ToField2(gene_id_list, "id", "symbol"))
    },

    symbolToId = function(gene_symbol_list) {
      return(private$fromFiled1ToField2(gene_symbol_list, "symbol", "id"))
    },

    idToEnsemblId = function(gene_id_list) {
      return(private$fromFiled1ToField2(gene_id_list, "id", "ensemblId"))
    },

    ensemblIdToId = function(gene_ensembl_id_list) {
      return(private$fromFiled1ToField2(gene_ensembl_id_list, "ensemblId", "id"))
    },

    ensemblIdToSymbol = function(gene_ensembl_id_list) {
      return(private$fromFiled1ToField2(gene_ensembl_id_list, "ensemblId", "symbol"))
    },

    symbolToEnsemblId = function(gene_symbol_list) {
      return(private$fromFiled1ToField2(gene_symbol_list, "symbol", "ensemblId"))
    }
  ),
  private = list(
    gene_set = NA,
    fromFiled1ToField2 = function(gene_list, from_field_name, to_field_name) {
      select_genes <- subset(private$gene_set, private$gene_set[[from_field_name]] %in% gene_list)
      select_genes <- select_genes[match(gene_list, select_genes[[from_field_name]]),]
      return(as.character(select_genes[[to_field_name]]))
    }
  )
)
