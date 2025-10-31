source("modules/disease_signature/GeneFilterByProteinCoding.R")


geneFilterByProteinCoding <- GeneFilterByProteinCoding$new()

gene_expressions_1 <- read.table("output/disease_signature/ipf/ipf_bin_chen_most_significant_genes.csv", sep = ";", header = T)
rownames(gene_expressions_1) <- gene_expressions_1$gene
result_1 <- geneFilterByProteinCoding$filter(gene_expressions_1)

gene_expressions_2 <- read.table("output/disease_signature/ipf/ipf_150_most_significant_genes.csv", sep = ";", header = T)
rownames(gene_expressions_2) <- gene_expressions_2$gene
result_2 <- geneFilterByProteinCoding$filter(gene_expressions_2)

gene_expressions_5 <- read.table("output/disease_signature/ipf/ipf_150_most_significant_landmark_genes.csv", sep = ";", header = T)
rownames(gene_expressions_5) <- gene_expressions_5$gene
result_5 <- geneFilterByProteinCoding$filter(gene_expressions_5)

gene_expressions_3 <- read.table("output/disease_signature/als_1/als_1_150_most_significant_genes.csv", sep = ";", header = T)
rownames(gene_expressions_3) <- gene_expressions_3$gene
result_3 <- geneFilterByProteinCoding$filter(gene_expressions_3)

gene_expressions_4 <- read.table("output/disease_signature/als_1/als_1_150_most_significant_landmark_genes.csv", sep = ";", header = T)
rownames(gene_expressions_4) <- gene_expressions_4$gene
result_4 <- geneFilterByProteinCoding$filter(gene_expressions_4)