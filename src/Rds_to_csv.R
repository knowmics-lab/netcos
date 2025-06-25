# Created 18/01/2025 to convert LINCS/metanalysis into .csv in local
# used to test my python versoin of the conn_score calculator, and to test it against
# original DEG data, in which the metanalysis were only Rds, which cannot be read by python.
# (MIGHT NEED TO CREATE A .pkl version too if it's too slow with csvs..)
setwd("~/unict_2024-25/drug_repurposing")
library(tools)
gene_wise_metanalysis_path='tsr/output/drug_signature/LINCS/metanalysis/'



for (gene in list.files(gene_wise_metanalysis_path))
  {
    gene_file=readRDS(paste(gene_wise_metanalysis_path,gene, sep=''))
    filename=paste(gene_wise_metanalysis_path,tools::file_path_sans_ext(gene),'.csv','',sep='')
    write.table(gene_file, filename, row.names=FALSE, col.names=TRUE, sep=';', quote=FALSE, dec = ",")
}

