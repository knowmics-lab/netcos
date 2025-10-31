DichotomicVoomDGE <- R6Class(
  "DichotomicVoomDGE",
  public = list(
    initialize = function(geneFilter = NA) {
      private$dicotomicRNASeqMapper <- DicotomicRNASeqMapper$new(geneFilter)
      private$voomDGEMapper <- VoomDGEMapper$new()
      private$geneFilterByProteinCoding <- GeneFilterByProteinCoding$new()
    },
    compute = function(rna_seq, sample_01_map, test_sample_name, filter_by_protein_coding = F) {
      gene_experiments_data <- private$dicotomicRNASeqMapper$load(rna_seq, sample_01_map, test_sample_name)
      if (filter_by_protein_coding) {
        gene_experiments_data$gene_expressions <- private$geneFilterByProteinCoding$filterById(gene_experiments_data$gene_expressions)
      }
      samples_metadata <- data.frame(
        sample_types = gene_experiments_data$sample_types,
        sample = colnames(gene_experiments_data$gene_expressions)
      )
      dge <- DGEList(gene_experiments_data$gene_expressions, remove.zeros = TRUE)
      dge <- calcNormFactors(dge, method = 'upperquartile')
      design <- model.matrix(~sample_types, data = samples_metadata)
      voom_data <- voom(dge, design, plot = F)
      fit_voom <- lmFit(voom_data, design)
      eBayes_fit_voom <- eBayes(fit_voom)
      differential_expression <- topTable(eBayes_fit_voom, coef = 2, number = 10^6)
      differential_expression$std.error <- differential_expression$logFC / differential_expression$t
      return(private$voomDGEMapper$map(differential_expression))
    }
  ),
  private = list(
    dicotomicRNASeqMapper = NA,
    geneFilterByProteinCoding = NA,
    voomDGEMapper = NA
  )
)
