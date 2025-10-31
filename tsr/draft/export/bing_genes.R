source("modules/config.R")
source("draft/export/load_LINCS_gene_info_list.R")
gene_info_list <- load_LINCS_gene_info_list(drug_signature_cfg$gene_info_filename)
gene_info_list <- subset(gene_info_list, gene_info_list$pr_is_bing == 1)
gene_info_list<-gene_info_list[,c("pr_gene_id","pr_gene_symbol","pr_is_lm","pr_is_bing")]
gene_info_list<-gene_info_list[order(gene_info_list$pr_gene_id),]
write.table(gene_info_list,"data/GSE92742-LINCS/bing_genes.csv",sep=";",row.names = F)