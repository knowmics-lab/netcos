ipf_chen_genes <- read.table("draft/data/ipf_chen_genes.csv", sep = ";", header = T)
ipf_150_most_significant_genes <- read.table("output/disease_signature/ipf/ipf_150_most_significant_genes.csv", sep = ";", header = T)

genes_to_compute<-subset(ipf_chen_genes, !ipf_chen_genes$gene %in% ipf_150_most_significant_genes$gene)
genes_to_compute$gene<-sort(genes_to_compute$gene)

write.table(genes_to_compute, "output/disease_signature/ipf/bin_chen_original_most_significant.csv", row.names = F, quote = F)