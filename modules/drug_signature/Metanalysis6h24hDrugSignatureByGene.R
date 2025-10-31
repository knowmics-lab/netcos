library(R6)
library(foreach)

source("modules/config.R")
source("modules/drug_signature/Metanalysis6h24hDrugSignature.R")
source("modules/TSRLogger.R")
source("modules/obj_is_na.R")

Metanalysis6h24hDrugSignatureByGene <- R6Class(
  "Metanalysis6h24hDrugSignatureByGene",
  public = list(
    initialize = function(drug_by_gene_signates_dir = NA, metanalysis_6h_24_dir = NA) {
      private$metanalysis6h24hDrugSignature <- Metanalysis6h24hDrugSignature$new()
      private$tsrLogger <- TSRLogger$new()
      if (obj_is_na(drug_by_gene_signates_dir))
        private$drug_by_gene_signates_dir <- drug_signature_cfg$genes_signature_full_path_dir
      else
        private$drug_by_gene_signates_dir <- drug_by_gene_signates_dir

      if (obj_is_na(metanalysis_6h_24_dir))
        private$metanalysis_6h_24_dir <- drug_signature_cfg$metanalysis_6h_24h_full_path_dir
      else
        private$metanalysis_6h_24_dir <- metanalysis_6h_24_dir
    },
    compute = function(drug_by_gene_signatures_base_filename, save_file = T) {
      startTime <- Sys.time()
      private$tsrLogger$log(sprintf("start meta analysis drug by gene signature 6h/24h: %s", drug_by_gene_signatures_base_filename))
      if (save_file) {
        filename <- paste0(private$metanalysis_6h_24_dir, drug_by_gene_signatures_base_filename, "_metanalysis.Rds")
        if (file.exists(filename)) {
          private$tsrLogger$log(sprintf("skipping gene %s", drug_by_gene_signatures_base_filename))
          return(NA)
        }
      }
      ### Combine 6h and 24h using meta-analysis: ###
      drug_by_gene_signatures_6h <- readRDS(paste0(private$drug_by_gene_signates_dir, drug_by_gene_signatures_base_filename, "_6h.Rds"))
      drug_by_gene_signatures_24h <- readRDS(paste0(private$drug_by_gene_signates_dir, drug_by_gene_signatures_base_filename, "_24h.Rds"))

      # calcola l'intersezione dei nomi dei farmaci tra gli insieme drugs_6h e drugs_24h.
      # merge fa la natural join
      drug_by_gene_signatures_6h_24h <- merge(x = drug_by_gene_signatures_6h, y = drug_by_gene_signatures_24h, by = c("gene", "drug"))
      # assegna i nomi alle colonne
      colnames(drug_by_gene_signatures_6h_24h)[3:12] <- c("DE_log2_FC_6h", "std.error_6h", "t.value_6h", "p.value_6h", "adj.p.value_6h")
      colnames(drug_by_gene_signatures_6h_24h)[8:12] <- c("DE_log2_FC_24h", "std.error_24h", "t.value_24h", "p.value_24h", "adj.p.value_24h")

      # in questo ciclo for viene fatta la meta analisi per calcolare l'effetto combinato di drugs_6h e drugs_24h
      total_drug_by_gene_signatures <- nrow(drug_by_gene_signatures_6h_24h)
      metalalysis_6h_24h <- data.frame(matrix(NA, nrow = total_drug_by_gene_signatures, ncol = 16))
      private$tsrLogger$log(sprintf("metanalysis to compute by gene %s: %s", drug_by_gene_signatures_base_filename, total_drug_by_gene_signatures))
      for (i in 1:total_drug_by_gene_signatures) {
        metalalysis_6h_24h[i,] <- private$metanalysis6h24hDrugSignature$compute(drug_by_gene_signatures_6h_24h[i,])
      }
      colnames(metalalysis_6h_24h) <- c(colnames(drug_by_gene_signatures_6h_24h), c("DE_log2_FC_6h_24h", "std.error_6h_24h", "t.value_6h_24h", "p.value_6h_24h"))
      # drug_by_gene_signatures_6h_24h_ma
      #    drug      gene   DE_log2_FC_6h    std.error_6h   t.value_6h      DE_log2_FC_24h    std.error_24h        t.value_24h      DE_log2_FC_6h_24h  std.error_6h_24h    t.value_6h_24h       p.value_6h_24h
      # 1  ABT-737   AARS   0.4461631309544  0.14302313904  0.155102125097  0.15510212509716  0.15479839964925124  1.0019620709813  0.3066570115074    0.1454057567166199  0.21089743517179400  0.0349467955686023
      # 2  ABT-751   AARS   0.0103021782160  0.05283095800  0.195002676573  -0.0534011868235  0.09920129528376088  -0.538311386668  -0.003773423983    0.0466304544571578  -0.0809218787845287  0.9355040802701673
      if (save_file) {
        saveRDS(metalalysis_6h_24h, filename)
      }
      totalTime <- Sys.time() - startTime
      private$tsrLogger$log(sprintf("end meta analysis drug by gene signature 6h/24h: %s, computation time: %s %s", drug_by_gene_signatures_base_filename, totalTime, attr(totalTime, "units")))
      return(metalalysis_6h_24h)
    }
  ),
  private = list(
    metanalysis6h24hDrugSignature = NA,
    tsrLogger = NA,
    drug_by_gene_signates_dir = NA,
    metanalysis_6h_24_dir = NA
  )
)