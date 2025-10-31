library("readr")

gene_symbols <- read_delim("data/LINCS-GSE92742/GSE92742_Broad_LINCS_gene_info.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

bing_gene_symbols <- subset(gene_symbols, gene_symbols$pr_is_bing == 1)

bing_gene_symbols <- sort(bing_gene_symbols$pr_gene_symbol)

write.table(bing_gene_symbols, "data/LINCS-GSE92742/bing_gene_symbols.csv", row.names = F, col.names = "gene", quote = F)
