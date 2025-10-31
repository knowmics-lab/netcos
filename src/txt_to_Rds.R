# 
# converts all fles in target folder from .txt to .Rds
#

setwd("~/unict_2024-25/drug_repurposing")
# choose directory and extension
# this will create Rds files for all files in given directory, for given extension
#target_directory='tsr/output/drug_signature/LINCS_lorenzo/metanalysis_mith3_gene_wise/'
#pattern='*.txt'
target_directory='tsr/output/disease_signature/als_NYGC_lorenzo/'
pattern='*.csv'
#

files_list= list.files(path=target_directory, pattern=pattern, full.names=FALSE, recursive=FALSE)

# save file
for (file_name in files_list){
data=read.table(paste(target_directory, file_name, sep=''), header = TRUE, sep = "\t", dec = ".")
rds_filename= substr(file_name,1,nchar(file_name)-4)
if (! file.exists(paste(target_directory,rds_filename,'.Rds',sep='')))
{
  print(file_name)
  saveRDS(data, paste(target_directory,rds_filename,'.Rds',sep=''))
}
}

