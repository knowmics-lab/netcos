gene_symbols <- read.table("data/LINCS-GSE92742/bing_gene_symbols.csv", header = T)$gene

gene_symbols<-sort(gene_symbols)

size <- length(gene_symbols) %/% 4

gene_split_1 <- gene_symbols[1:size]
write.table(gene_split_1, "data/LINCS-GSE92742/bing_gene_symbols_split_1.csv", row.names = F, col.names = "gene", quote = F)
gene_split_2 <- gene_symbols[(size + 1):(2 * size)]
write.table(gene_split_2, "data/LINCS-GSE92742/bing_gene_symbols_split_2.csv", row.names = F, col.names = "gene", quote = F)
gene_split_3 <- gene_symbols[(2 * size + 1):(3 * size)]
write.table(gene_split_3, "data/LINCS-GSE92742/bing_gene_symbols_split_3.csv", row.names = F, col.names = "gene", quote = F)
gene_split_4 <- gene_symbols[(3 * size + 1):length(gene_symbols)]
write.table(gene_split_4, "data/LINCS-GSE92742/bing_gene_symbols_split_4.csv", row.names = F, col.names = "gene", quote = F)
