library(readr)
colon_tumor_signature <- read_delim("data/colon_tumor/colon_tumor_human_signature.tsv", delim = "\t")
write.table(colon_tumor_signature, "data/colon_tumor/colon_tumor_human_signature.csv", sep = ";", dec = ",", row.names = F)
lincs_gene_symbols <- read.table("data/LINCS-GSE92742/bing_genes.csv", sep = ";", header = T)
colnames(lincs_gene_symbols)[2] <- "gene"

common_genes <- merge(x = colon_tumor_signature, y = lincs_gene_symbols, by.x = "Gene Name", by.y = "gene")
common_genes <- common_genes[order(common_genes$adj.p, decreasing = F),]

common_genes$`Gene Id` <- NULL

colnames(common_genes) <- c("gene", "DE_log2_FC", "p.value", "adj.p.value", "gene_id", "is_landmark", "is_bing")

common_genes <- common_genes[, c("gene", "DE_log2_FC", "DE_log2_FC", "p.value", "adj.p.value", "gene_id", "is_landmark", "is_bing")]
colnames(common_genes)[3] <- "t.value"
#significant_genes <- subset(common_genes, adj.p.value < 0.001)
#significant_genes <- significant_genes[order(abs(significant_genes$DE_log2_FC), decreasing = T),]
#write.table(significant_genes, "data/colon_tumor/colon_tumor_significant_genes.csv", row.names = F, quote = F, sep = ";", dec = ",")

#significant_genes$abslog<-abs(significant_genes$DE_log2_FC)
#colon_tumor_150_most_significant_genes <- significant_genes[1:150,]
#saveRDS(colon_tumor_150_most_significant_genes, "data/colon_tumor/colon_tumor_150_most_significant_genes.Rds")
#write.table(colon_tumor_150_most_significant_genes, "data/colon_tumor/colon_tumor_150_most_significant_genes.csv", row.names = F, quote = F, sep = ";", dec = ",")

landmark_genes <- subset(common_genes, common_genes$is_landmark == 1)

colon_tumor_150_most_significant_landmark_genes <- landmark_genes[1:151,]
colon_tumor_150_most_significant_landmark_genes <- colon_tumor_150_most_significant_landmark_genes[!duplicated(colon_tumor_150_most_significant_landmark_genes$gene),]
saveRDS(colon_tumor_150_most_significant_landmark_genes, "data/colon_tumor/colon_tumor_150_most_significant_landmark_genes.Rds")
write.table(colon_tumor_150_most_significant_landmark_genes, "data/colon_tumor/colon_tumor_150_most_significant_landmark_genes.csv", row.names = F, quote = F, sep = ";", dec = ",")

print("data written")
