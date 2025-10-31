source("modules/config.R")
ipf_signature <- readRDS("output/disease_signature/ipf/ipf_signature.Rds")
lincs_gene_symbols <- read.table("data/LINCS-GSE92742/bing_genes.csv", sep = ";", header = T)
colnames(lincs_gene_symbols)[2] <- "gene"

common_genes <- merge(x = ipf_signature, y = lincs_gene_symbols, by.x = "gene", by.y = "gene")
common_genes$abs_t.value <- abs(common_genes$t.value)
common_genes <- common_genes[order(common_genes$abs_t.value, decreasing = T),]
colnames(common_genes) <- c("gene", "DE_log2_FC", "std.error", "t.value", "p.value", "adj.p.value", "gene_id", "is_landmark", "is_bing", "abs_t.value")
common_genes <- common_genes[, c("gene", "DE_log2_FC", "std.error", "t.value", "abs_t.value", "p.value", "adj.p.value", "gene_id", "is_landmark", "is_bing")]

ipf_bin_chen_most_significant <- subset(common_genes, common_genes$adj.p.value <= 0.001 & abs(common_genes$DE_log2_FC) > 1.5)
ipf_bin_chen_most_significant <- ipf_bin_chen_most_significant[order(ipf_bin_chen_most_significant$gene),]
saveRDS(ipf_bin_chen_most_significant, gsub(".csv", ".Rds", ipf_disease_signature_cfg$ipf_bin_chen_most_significant_genes_filename))
write.table(ipf_bin_chen_most_significant, ipf_disease_signature_cfg$ipf_bin_chen_most_significant_genes_filename, row.names = F, quote = F, sep = ";", dec = ",")

ipf_150_most_significant_genes <- common_genes[1:150,]
saveRDS(ipf_150_most_significant_genes, gsub(".csv", ".Rds", ipf_disease_signature_cfg$ipf_150_most_significant_genes_filename))
write.table(ipf_150_most_significant_genes, ipf_disease_signature_cfg$ipf_150_most_significant_genes_filename, row.names = F, quote = F, sep = ";", dec = ",")

landmark_genes <- subset(common_genes, common_genes$is_landmark == 1)

ipf_150_most_significant_landmark_genes <- landmark_genes[1:150,]
saveRDS(ipf_150_most_significant_landmark_genes, gsub(".csv", ".Rds", ipf_disease_signature_cfg$ipf_150_most_significant_landmark_genes_filename))
write.table(ipf_150_most_significant_landmark_genes, ipf_disease_signature_cfg$ipf_150_most_significant_landmark_genes_filename, row.names = F, quote = F, sep = ";", dec = ",")

print("data written")
