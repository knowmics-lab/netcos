source("modules/config.R")
source("modules/drug_signature/setuper/LINCSExperimentMetaDataSetuper.R")
source("modules/drug_signature/loader/LINCSExperimentDataRowsLoader.R")
source("modules/drug_signature/exporter/LINCSExporterByGeneList.R")

LINCSExporter <- R6Class(
  "LINCSExporter",
  public = list(
    initialize = function(lincs_splitted_level3_dir = NA) {
      private$lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)
      lincsExperimentDataRowsLoader <- LINCSExperimentDataRowsLoader$new(drug_signature_cfg$experiments_data_filename)
      private$lincsExporterByGeneList <- LINCSExporterByGeneList$new(lincsExperimentDataRowsLoader, lincs_splitted_level3_dir)
    },
    export = function(filename, perturbation_times = drug_signature_cfg$perturbation_times) {
      genes <- read.table(filename, sep = ";", header = T)
      genes <- genes[order(genes$gene_id),]
      experiments_env <- private$lincsExperimentMetaDataSetuper$setup(perturbation_times)
      block_size <- 100
      gene_index <- 1
      tot_genes <- length(genes$gene_id)
      print(sprintf("export %s genes", tot_genes))
      while (gene_index <= tot_genes) {
        ext_sup <- min(gene_index + block_size - 1, tot_genes)
        private$lincsExporterByGeneList$export(genes, experiments_env, gene_index, ext_sup)
        gene_index <- gene_index + block_size
      }
    }
  ),
  private = list(
    lincsExperimentMetaDataSetuper = NA,
    lincsExporterByGeneList = NA
  )
)