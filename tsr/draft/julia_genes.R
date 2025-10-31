lincs_genes <- read.table("data/LINCS-GSE92742/bing_genes.csv", sep = ";", header = T)
colnames(lincs_genes)[2] <- "gene"
colnames(lincs_genes)[1] <- "id"
ipf_150_most_significant_genes <- read.table("output/disease_signature/ipf/ipf_150_most_significant_genes.csv", sep = ";", header = T)

lincs_genes <- subset(lincs_genes, !lincs_genes$gene %in% ipf_150_most_significant_genes$gene)

genes <- NA
filenames <- NA
for (i in 1:dim(lincs_genes)[1]) {
  genes[i * 2 - 1] <- lincs_genes$gene[i]
  filenames[i * 2 - 1] <- paste0(lincs_genes$gene[i], "_", lincs_genes$id[i], "_6h")
  genes[i * 2] <- lincs_genes$gene[i]
  filenames[i * 2] <- paste0(lincs_genes$gene[i], "_", lincs_genes$id[i], "_24h")
}

export_genes <- data.frame(
  gene = genes,
  filename = filenames
)
write.table(export_genes, "data/LINCS-GSE92742/julia_genes_to_compute.csv", sep = ";", row.names = F, quote = F)