source("modules/drug_signature/converter/GeneSymbolToIdConverter.R")
source("modules/drug_signature/loader/GeneInfoListLoader.R")

geneSymbolToIdConverter <- GeneSymbolToIdConverter$new(GeneInfoListLoader$new())

colon_tumor_human_genes <- read.table("draft/data/colon_tumor_human.csv", header = T, sep = ";", dec=",")
for (i in 1:dim(colon_tumor_human_genes)[1]) {
  colon_tumor_human_genes$gene[i] <- toupper(colon_tumor_human_genes$gene[i])
}

colon_tumor_mouse_genes <- read.table("draft/data/colon_tumor_mouse.csv", header = T, sep = ";", dec=",")
for (i in 1:dim(colon_tumor_mouse_genes)[1]) {
  colon_tumor_mouse_genes$gene[i] <- toupper(colon_tumor_mouse_genes$gene[i])
}

bing_genes <- read.table("data/LINCS-GSE92742/bing_gene_symbols.csv", sep = ";", header = T)

common_colon_tumor_human_genes <- subset(colon_tumor_human_genes, gene %in% bing_genes$gene)
human_gene_ids <- geneSymbolToIdConverter$convert(common_colon_tumor_human_genes$gene)
colon_tumor_human_signature <- data.frame(
  gene = common_colon_tumor_human_genes$gene,
  DE_log2_FC = common_colon_tumor_human_genes$DE_log2_FC,
  std.error = NA,
  t.value = common_colon_tumor_human_genes$DE_log2_FC,
  p.value = NA,
  adj.p.value = NA,
  gene_id = human_gene_ids
)
write.table(colon_tumor_human_signature, "draft/data/colon_tumor_human_signature.csv", row.names = F, sep = ";", dec = ",")
saveRDS(colon_tumor_human_signature,"draft/data/colon_tumor_human_signature.Rds")

common_colon_tumor_mouse_genes <- subset(colon_tumor_mouse_genes, gene %in% bing_genes$gene)
mouse_gene_ids <- geneSymbolToIdConverter$convert(common_colon_tumor_mouse_genes$gene)
colon_tumor_mouse_signature <- data.frame(
  gene = common_colon_tumor_mouse_genes$gene,
  DE_log2_FC = common_colon_tumor_mouse_genes$DE_log2_FC,
  std.error = NA,
  t.value = common_colon_tumor_mouse_genes$DE_log2_FC,
  p.value = NA,
  adj.p.value = NA,
  gene_id = mouse_gene_ids
)
write.table(colon_tumor_mouse_signature, "draft/data/colon_tumor_mouse_signature.csv", row.names = F, sep = ";", dec = ",")
saveRDS(colon_tumor_mouse_signature,"draft/data/colon_tumor_mouse_signature.Rds")

colon_tumor_signature <- colon_tumor_human_signature
colon_tumor_signature <- rbind(colon_tumor_human_signature, colon_tumor_mouse_signature)
genes<-unique(colon_tumor_signature$gene)

write.table(genes, "draft/data/colon_tumor_genes.csv", row.names = F, sep = ";", dec = ",", col.names = "gene")



print("data written successfully")