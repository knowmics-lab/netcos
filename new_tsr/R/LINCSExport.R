
LINCSExport <- R6Class(
  "LINCSExport",
  public = list(
    initialize = function(experiments_data_filename, lincsMetadataSetuper = NA) {
      if (obj_is_na(lincsMetadataSetuper)){
        private$lincsMetadataSetuper <- LINCSMetadataSetuper$new()
      }else {
        if (!"LINCSMetadataSetuperAbstract" %in% class(lincsMetadataSetuper))
          stop("the lincsExperimentMetaDataSetuper instance must be of type LINCSMetadataSetuperAbstract")
        private$lincsMetadataSetuper <- lincsMetadataSetuper
      }
      lincsSGCTXDataRowLoader <- LINCSGCTXDataRowLoader$new(experiments_data_filename)
      private$lincsExportByGeneList <- LINCSExportByGeneList$new(lincsSGCTXDataRowLoader)
    },
    export = function(gene_ids,
                      perturbation_times,
                      output_dir,
                      block_size = 200,
                      drugs_filter = NA
    ) {
      experiments_meta_data <- private$lincsMetadataSetuper$setup(perturbation_times, drugs_filter)
      gene_index <- 1
      tot_genes <- length(gene_ids)
      while (gene_index <= tot_genes) {
        ext_sup <- min(gene_index + block_size - 1, tot_genes)
        private$lincsExportByGeneList$export(gene_ids, experiments_meta_data, gene_index, ext_sup, output_dir)
        gene_index <- gene_index + block_size
      }
    }
  ),
  private = list(
    lincsMetadataSetuper = NA,
    lincsExportByGeneList = NA
  )
)
