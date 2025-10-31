tcga_csv <- read.table("draft/data/disease_signature_TCGA_for_test.csv",sep = ";",header = T)
saveRDS(tcga_csv,"draft/data/disease_signature_TCGA_for_test.Rds")
print("ok")
