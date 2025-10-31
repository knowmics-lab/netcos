library(R6)
source("modules/config.R")
source("modules/TSRLogger.R")
MostSignificantGenes <- R6Class(
  "MostSignificantGenes",
  public = list(
    compute = function(disease_name, signature_filename, bin_chen_most_significant_filename, a150_most_significant_genes_filename, a150_most_significant_landmark_genes_filename, bin_chen_min_FC = 1.5) {
      private$tsrLogger$log(sprintf("start %s most significant genes computation", disease_name))

      signature <- readRDS(signature_filename)
      lincs_gene_symbols <- read.table(private$bing_genes_filename, sep = ";", header = T)
      colnames(lincs_gene_symbols)[2] <- "gene"

      common_genes <- merge(x = signature, y = lincs_gene_symbols, by.x = "gene", by.y = "gene")
      common_genes$abs_t.value <- abs(common_genes$t.value)
      common_genes <- common_genes[order(common_genes$abs_t.value, decreasing = T),]
      colnames(common_genes) <- c("gene", "DE_log2_FC", "std.error", "t.value", "p.value", "adj.p.value", "gene_id", "is_landmark", "is_bing", "abs_t.value")
      common_genes <- common_genes[, c("gene", "DE_log2_FC", "std.error", "t.value", "abs_t.value", "p.value", "adj.p.value", "gene_id", "is_landmark", "is_bing")]

      bin_chen_most_significant <- subset(common_genes, common_genes$adj.p.value <= 0.001 & abs(common_genes$DE_log2_FC) > bin_chen_min_FC)
      bin_chen_most_significant <- bin_chen_most_significant[order(bin_chen_most_significant$gene),]
      saveRDS(bin_chen_most_significant, gsub(".csv", ".Rds", bin_chen_most_significant_filename))
      write.table(bin_chen_most_significant, bin_chen_most_significant_filename, row.names = F, quote = F, sep = ";", dec = ",")

      a150_most_significant_genes <- common_genes[1:150,]
      saveRDS(a150_most_significant_genes, gsub(".csv", ".Rds", a150_most_significant_genes_filename))
      write.table(a150_most_significant_genes, a150_most_significant_genes_filename, row.names = F, quote = F, sep = ";", dec = ",")

      landmark_genes <- subset(common_genes, common_genes$is_landmark == 1)

      a150_most_significant_landmark_genes <- landmark_genes[1:150,]
      saveRDS(a150_most_significant_landmark_genes, gsub(".csv", ".Rds", a150_most_significant_landmark_genes_filename))
      write.table(a150_most_significant_landmark_genes, a150_most_significant_landmark_genes_filename, row.names = F, quote = F, sep = ";", dec = ",")

      private$tsrLogger$log(sprintf("end %s most signafucant genes computation", disease_name))

    }
  ),
  private = list(
    bing_genes_filename = drug_signature_cfg$bing_genes_filename,
    tsrLogger = TSRLogger$new()
  )
)