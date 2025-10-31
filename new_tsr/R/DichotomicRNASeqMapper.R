DicotomicRNASeqMapper <- R6Class(
  "DicotomicRNASeqMapper",
  public = list(
    initialize = function(geneFilter = NA) {
      if (!obj_is_na(geneFilter)) {
        if (!"GeneFilterAbstract" %in% class(geneFilter))
          stop("the geneFilter instance must by of type GeneFilterAbstract")
        private$geneFilter <- geneFilter
      }else {
        private$geneFilter <- LowCountsGeneFilter$new()
      }
      private$rnaSeqSampleMapBuilder <- RnaSeqSampleMapBuilder$new()

    },
    load = function(rna_seq, sample_01_map, test_sample_name) {
      # rna_seq = matrix, rownames = gene_id, cols=experiment rna seq (read count)
      # Example, GSM2433098, GSM2433099, ... = colnames;  100287102, 653635, ...= rownames
      # colnames are not used so they are optional
      #     	    GSM2433098	GSM2433099	GSM2433100	GSM2433101	GSM2433102
      # 100287102	7	        4	        6	        7	        12
      # 653635	    817			480			513			497			1055
      # 102466751	30			18			14			20			24
      # 107985730	1			0			0			0			2
      rnaSeqSampleMap <- private$rnaSeqSampleMapBuilder$build(sample_01_map, test_sample_name)
      rna_seq <- rna_seq[, rnaSeqSampleMap$sample_positions]
      rna_seq <- private$geneFilter$filter(rna_seq, rnaSeqSampleMap$samples)
      return(
        list(
          gene_expressions = rna_seq,
          sample_types = rnaSeqSampleMap$samples
        )
      )
    }
  ),
  private = list(
    rnaSeqSampleMapBuilder = NA,
    geneFilter = NA
  )
)