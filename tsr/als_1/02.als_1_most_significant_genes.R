source("modules/config.R")
#als_1_signature <- readRDS("output/disease_signature/als_1/als_signature.Rds")
als_1_signature <- read.table("output/disease_signature/als_1/als_signature_fixed.csv", sep = ";", header = T, dec = ",")
lincs_gene_symbols <- read.table("data/LINCS-GSE92742/bing_genes.csv", sep = ";", header = T)
colnames(lincs_gene_symbols)[2] <- "gene"

common_genes <- merge(x = als_1_signature, y = lincs_gene_symbols, by.x = "gene", by.y = "gene")
common_genes$abs_t.value <- abs(common_genes$t.value)
common_genes <- common_genes[order(common_genes$abs_t.value, decreasing = T),]
colnames(common_genes) <- c("gene", "DE_log2_FC", "std.error", "t.value", "p.value", "adj.p.value", "gene_id", "is_landmark", "is_bing", "abs_t.value")
common_genes <- common_genes[, c("gene", "DE_log2_FC", "std.error", "t.value", "abs_t.value", "p.value", "adj.p.value", "gene_id", "is_landmark", "is_bing")]

als_1_bin_chen_most_significant <- subset(common_genes, common_genes$adj.p.value <= 0.001 & abs(common_genes$DE_log2_FC) > 1.5)
als_1_bin_chen_most_significant <- als_1_bin_chen_most_significant[order(als_1_bin_chen_most_significant$gene),]
write.table(als_1_bin_chen_most_significant, als_1_disease_signature_cfg$als_1_bin_chen_most_significant_genes_filename, row.names = F, quote = F, sep = ";", dec = ",")

als_1_150_most_significant_genes <- common_genes[1:150,]
write.table(als_1_150_most_significant_genes, als_1_disease_signature_cfg$als_1_150_most_significant_genes_filename, row.names = F, quote = F, sep = ";", dec = ",")

landmark_genes <- subset(common_genes, common_genes$is_landmark == 1)

als_1_150_most_significant_landmark_genes <- landmark_genes[1:150,]
write.table(als_1_150_most_significant_landmark_genes, als_1_disease_signature_cfg$als_1_150_most_significant_landmark_genes_filename, row.names = F, quote = F, sep = ";", dec = ",")

print("data written")
