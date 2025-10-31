setwd("G:/My Drive/unict 2024-25/drug_repurposing")
source("modules/config.R")
source("modules/drug_signature/setuper/LINCSExperimentMetaDataSetuper.R")
source("modules/LINCS_export/LINCSExperimentDataRowsLoader.R")
source("modules/LINCS_export/LINCSExportByGeneList.R")

lincsExperimentMetaDataSetuper <- LINCSExperimentMetaDataSetuper$new(drug_signature_cfg$experiments_meta_data_filename)
lincsExperimentDataRowsLoader <- LINCSExperimentDataRowsLoader$new(drug_signature_cfg$experiments_data_filename)
lincsExportByGeneList <- LINCSExportByGeneList$new(lincsExperimentDataRowsLoader, drug_signature_cfg$lincs_splitted_level3_dir)

bing_genes <- read.table("data/LINCS-GSE92742/bing_genes.csv", sep = ";", header = T)
experiments_env <- lincsExperimentMetaDataSetuper$setup(c("6", "24"))
block_size <- 400
gene_index <- 1
tot_genes <- length(bing_genes$pr_gene_id) # 10174
print(sprintf("computation from index: %s to index: %s", gene_index, tot_genes))
while (gene_index <= tot_genes) {
  ext_sup <- min(gene_index + block_size - 1, tot_genes)
  lincsExportByGeneList$export(bing_genes, experiments_env, gene_index, ext_sup)
  gene_index <- gene_index + block_size
}
