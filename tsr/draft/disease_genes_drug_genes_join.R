drug_genes <-
  read_delim(
    "data/LINCS/GSE92742_Broad_LINCS_gene_info.txt",
    "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )

disease_genes <- readRDS("output/disease_signature/als/als_signature.Rds")

common_genes <- merge(x = disease_genes, y = drug_genes, by.x = "gene", by.y = "pr_gene_symbol")
write.table(common_genes[, "gene"], "data/common_genes.csv", sep = ";", row.names = F)
drug_landmark_genes <- subset(drug_genes, pr_is_lm == 1)

common_genes_landmark <- merge(x = disease_genes, y = drug_landmark_genes, by.x = "gene", by.y = "pr_gene_symbol")
write.table(common_genes_landmark[, "gene"], "data/common_genes_landmark.csv", sep = ";", row.names = F)

print("end computation")

