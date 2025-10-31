library(GEOquery)
library(R6)
source("modules/disease_signature/als_1/ALSInsignificantGeneFilter.R")
source("modules/disease_signature/als_1/ALSGeneExpressionFormatter.R")
source("modules/disease_signature/SamplesGroups.R")
source("modules/config.R")

ALS_GSE3307_GPL97DataPrepare <- R6Class(
  "ALS_GSE3307_GPL97DataPrepare",
  public = list(
    prepare = function(gse_3307_GPL97_filename, gse_3307_GPL97_temp_annot_dir, samples_groups_map) {
      gse_3307_gpl97 <- getGEO(filename = gse_3307_GPL97_filename, destdir = gse_3307_GPL97_temp_annot_dir, GSEMatrix = TRUE, AnnotGPL = TRUE)
      gene_expressions <- exprs(gse_3307_gpl97)
      als_samples_groups <- private$samplesGroups$groups(samples_groups_map, private$disease_name)
      gene_expressions <- gene_expressions[, als_samples_groups$sample_groups_idxes]
      gene_expressions <- private$alsGeneExpressionFormatter$format(gene_expressions, gse_3307_gpl97)
      gene_expressions <- private$insignificantGeneFilter$filter(gene_expressions)
      return(
        list(
          gene_expressions = gene_expressions,
          disease_control_groups = als_samples_groups$disease_control_groups
        )
      )
    }
  ),
  private = list(
    samplesGroups = SamplesGroups$new(),
    alsGeneExpressionFormatter = ALSGeneExpressionFormatter$new(),
    insignificantGeneFilter = ALSInsignificantGeneFilter$new(),
    disease_name = als_1_disease_signature_cfg$name
  )
)