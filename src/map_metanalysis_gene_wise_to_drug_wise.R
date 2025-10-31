setwd("/home/NetCos/drug_repurposing")
gene_wise_metanalysis_path='tsr/output/drug_signature/LINCS/metanalysis/'
drug_wise_metanalysis_path='tsr/output/drug_signature/LINCS/metanalysis_drug_wise/'
mithril_input_path='MITHRIL/input/durg_signature'
example_metanalysis=readRDS(paste(gene_wise_metanalysis_path,'BRCA1_672_metanalysis.Rds',sep=''))
example_drug_wise_metanalysis=readRDS(paste(drug_wise_metanalysis_path,'amitriptyline_metanalysis.Rds',sep=''))

drug_list=example_metanalysis$drug #  all files have same amount of drugs. should be around 3222
print(paste(length(drug_list),'drugs')) #should be 3222
print(names(example_metanalysis))

N=length(list.files(gene_wise_metanalysis_path))
print(paste(N,'genes'))

args=commandArgs(TRUE) # pass 2 command arguments corresponding to initial and final index of slice to slice drug_list

for (drug in drug_list[args[1]:args[2]])
{print(drug) 
  
  # Initialize empty dataframes for all timepoints 
  df_6=data.frame()
  df_24=data.frame()
  df_6_24=data.frame()
  df_drug_wise_metanalysis=data.frame()

  # Extract data for given timepoint from every gene file: 

  for (gene in list.files(gene_wise_metanalysis_path))
  {
    gene_file=readRDS(paste(gene_wise_metanalysis_path,gene, sep=''))
    gene_file$d
    # Select row corresponding to drug
    drug_row=gene_file[gene_file$drug==drug,]
    # Append it to drug_wise metanalysis dataframe
    df_drug_wise_metanalysis=rbind(df_drug_wise_metanalysis, drug_row) # the data.frame will take the column names of gene_file automatically
    
    # Append specific columns to mithril input dataframes
    df_6=rbind(df_6, drug_row[c("gene","DE_log2_FC_6h")])
    df_24=rbind(df_24, drug_row[c("gene","DE_log2_FC_24h")])
    df_6_24=rbind(df_6_24, drug_row[c("gene","DE_log2_FC_6h_24h")])
    
  
  }
  
  
  # Save drug-wise metanalysis rds file:
  print('saving drug-wise metanalysis Rds file')
  filename=paste(drug_wise_metanalysis_path,drug,'_metanalysis.Rds',sep='')
  saveRDS(df_drug_wise_metanalysis, filename)
  
  # Save MIThRIL input
  print('saving mithril inputs format for 6h, 24h and 6h_24h')

  filename=paste(mithril_input_path,drug,'_6h.mi',sep='')
  write.table(df_6, filename, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
  filename=paste(mithril_input_path,drug,'_24h.mi',sep='')
  write.table(df_24, filename, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
  filename=paste(mithril_input_path,drug,'_6h_24h.mi',sep='')
  write.table(df_6_24, filename, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

  
  }

  
