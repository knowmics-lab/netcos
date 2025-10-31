config <- list()
config$package_name <- "DrugRepurposing"
config$Illumina_HiSeq_annot_filename <- "extdata/Illumina_HiSeq/Human.GRCh38.p13.annot.tsv.gz"
config$protein_coding_gene_filename <- "extdata/protein_coding_gene.Rds"


config$gene_set_base_path <- "extdata/gene_sets/"
config$gene_set_filename_suffix <- "_gene_set.Rds"
config$gene_set_Human_GRCh38.p13 <- "Human_GRCh38.p13"
config$gene_set_LINCS <- "LINCS"
config$gene_set_hgnc <- "hgnc"
config$Illumina_HiSeq_annot_filename <- "extdata/Illumina_HiSeq/Human.GRCh38.p13.annot.tsv.gz"
config$protein_coding_genes_filename <- "extdata/protein_coding_genes.Rds"
config$LINCS_metadata_RDS_filename <- "extdata/GSE92742_Broad_LINCS_inst_info.Rds"
config$LINCSLMMFormula <- gene_expression ~ pert_iname + (1 | cell_id) + (1 | rna_plate)

