tcga <- readRDS("draft/data/disease_signature_TCGA.Rds")
meta <- readRDS("draft/data/gene_meta.Rds")
tcga_meta <- merge(x = tcga, y = meta, by.x = "gene", by.y = "Gene stable ID")
tcga_meta$gene<-conv$"Gene name"

new_tcga<-tcga_meta[,c("gene","Tumor", "log2_FC_estimate", "log2_FC_estimate_SE",
"t","P.Value", "adj.P.Val", "t.value" )]
saveRDS(new_tcga,"draft/data/disease_signature_TCGA2.Rds")
write.table(new_tcga,"draft/data/disease_signature_TCGA2.csv",sep = ";",row.names = F)
