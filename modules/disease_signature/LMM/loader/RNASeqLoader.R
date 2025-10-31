library(R6)
library(readr)

source("modules/TSRLogger.R")
source("modules/disease_signature/LMM/loader/RNASeqMetadataLoader.R")
source("modules/disease_signature/LMM/mapper/RNASeqMetadataMapper.R")
source("modules/disease_signature/LMM/loader/RNASeqDataLoader.R")
source("modules/disease_signature/Illumina_HiSeq_mapper/GeneExpressionFormatter.R")

RNASeqLoader <- R6Class(
  "RNASeqLoader",
  public = list(
    initialize = function(geneFilter) {
      if (!"GeneFilter" %in% class(geneFilter))
        stop("the geneFilter instance must by of type GeneFilter")
      private$geneFilter <- geneFilter
      private$rnaSeqMetadataLoader <- RNASeqMetadataLoader$new()
      private$rnaSeqMetadataMapper <- RNASeqMetadataMapper$new()
      private$rnaSeqDataLoader <- RNASeqDataLoader$new()
      private$geneExpressionFormatter <- GeneExpressionFormatter$new()
    },
    load = function(rna_seq_metadata_filename, rna_seq_data_filename, tissue_statuses_to_be_tested, tissue_statuses_map) {
      rna_seq_metadata <- private$rnaSeqMetadataLoader$load(rna_seq_metadata_filename)
      rna_seq_metadata <- private$rnaSeqMetadataMapper$map(rna_seq_metadata, tissue_statuses_to_be_tested, tissue_statuses_map)
      rna_seq_data <- private$rnaSeqDataLoader$load(rna_seq_data_filename, rna_seq_metadata)
      rna_seq_data <- private$geneFilter$filter(rna_seq_data, rna_seq_metadata$tissue_status)
      rna_seq_data <- private$geneExpressionFormatter$format(rna_seq_data)
      return(
        list(
          metadata = rna_seq_metadata,
          data = rna_seq_data
        )
      )
    }
  ),
  private = list(
    rnaSeqMetadataLoader = NA,
    rnaSeqMetadataMapper = NA,
    rnaSeqDataLoader = NA,
    geneFilter = NA,
    geneExpressionFormatter = NA
  )
)