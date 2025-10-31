library(R6)
library(lme4)
source("modules/drug_signature/mapper/abstract_class/LMMToDataFrameMapper.R")

LmerLMMToDataFrameMapper <- R6Class(
  "LmerLMMToDataFrameMapper",
  inherit = LMMToDataFrameMapper,
  public = list(
    map = function(LMM_output, gene) {
      # estrae i coefficienti dell'analisi "lmer". Esempio:
      #                               Estimate      Std. Error  t value
      # (Intercept)                   5.459338302   0.27436280  19.89824544
      # pert_inametyrphostin-AG-1296  -0.032810809  0.12512476  -0.26222474
      # pert_inametyrphostin-AG-1478  0.115920192   0.12512476  0.92643686
      gene_DE <- as.data.frame(coef(summary(LMM_output)))
      # crea la colonna drug dai nomi delle righe dei coefficienti dell'analisi "lmer"
      gene_DE$drug <- rownames(gene_DE)
      # elimina la prima riga contenente l'intercetta
      gene_DE <- gene_DE[2:nrow(gene_DE),]
      # elimina il prefisso pert_iname dal nome del farmaco nella colonna drug
      gene_DE$drug <- substr(gene_DE$drug, 11, nchar(gene_DE$drug))
      # crea una nuova colonna con il nome del gene. La colonna conterrà sempre lo stesso valore.
      gene_DE$gene <- gene
      # lascia solo le colonne contenenti i dati necessari
      gene_DE <- gene_DE[, c("drug", "gene", "Estimate", "Std. Error", "t value")]
      # rinomina le colonne
      colnames(gene_DE) <- c("drug", "gene", "DE_log2_FC", "std.error", "t.value")
      # elimina i nomi delle righe
      rownames(gene_DE) <- NULL
      # risultato finale:
      # drug                gene  Estimate      Std. Error  t value
      # tyrphostin-AG-1296  PAX8  -0.032810809  0.12512476  -0.26222474
      # tyrphostin-AG-1478  PAX8  0.115920192   0.12512476  0.92643686
      return(gene_DE)
    }
  )
)