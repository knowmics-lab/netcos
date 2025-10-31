lincs_gene <- read.table("data/LINCS-GSE92742/bing_genes.csv", sep = ";", header = T)
ipf_150_most_significant_genes <- read.table("output/disease_signature/ipf/ipf_150_most_significant_genes.csv", sep = ";", header = T)
lincs_gene_symbols <- lincs_gene[2]
colnames(lincs_gene_symbols) <- "gene"
lincs_gene_symbols <- subset(lincs_gene_symbols, !lincs_gene_symbols$gene %in% ipf_150_most_significant_genes$gene)
write.table(lincs_gene_symbols, "data/LINCS-GSE92742/genes_to_export.csv", row.names = F, quote = F)