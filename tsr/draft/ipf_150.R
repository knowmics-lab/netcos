ipf_gene_gymbols_signature <- readRDS("output/disease_signature/ipf/ipf_signature.Rds")
lincs_gene_symbols <- read.table("data/LINCS-GSE92742/bing_genes.csv", sep = ";", header = T)
colnames(lincs_gene_symbols)[2] <- "gene"

signatures <- merge(x = ipf_gene_gymbols_signature, y = lincs_gene_symbols, by.x = "gene", by.y = "gene")
signatures$abs_t.value <- abs(signatures$t.value)
signatures <- signatures[order(signatures$abs_t.value, decreasing = T),]

first_150_signatures <- signatures[1:150,]

signatures_landmark <- subset(signatures, signatures$pr_is_lm == 1)
signatures_landmark <- signatures_landmark[order(signatures_landmark$abs_t.value, decreasing = T),]
first_150_signatures_landmark <- signatures_landmark[1:150,]

first_genes <- rbind(first_150_signatures, first_150_signatures_landmark)
first_genes <- unique(first_genes)

first_l <- subset(first_150_signatures, first_150_signatures$pr_is_lm == 1)

gene_meta <- readRDS("data/gene_meta.Rds")

colnames(gene_meta)[7] <- "gene"

test <- merge(x = lincs_gene_symbols, y = gene_meta, by.x = "gene", by.y = "gene")
test2 <- subset(lincs_gene_symbols, !lincs_gene_symbols$gene %in% gene_meta$gene)
test3 <- subset(lincs_gene_symbols, !lincs_gene_symbols$gene %in% ipf_gene_gymbols_signature$gene)
test4 <- subset(ipf_gene_gymbols_signature, !ipf_gene_gymbols_signature$gene %in% lincs_gene_symbols$gene)