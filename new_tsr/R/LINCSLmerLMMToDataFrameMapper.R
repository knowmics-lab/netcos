
LINCSLmerLMMToDataFrameMapper <- R6Class(
  "LINCSLmerLMMToDataFrameMapper",
  inherit = LMMToDataFrameMapperAbstract,
  public = list(
    map = function(LMM_output, gene_id = NA) {
      # estrae i coefficienti dell'analisi "lmer". Esempio:
      #                               Estimate      Std. Error  t value
      # (Intercept)                   5.459338302   0.27436280  19.89824544
      # pert_inametyrphostin-AG-1296  -0.032810809  0.12512476  -0.26222474
      # pert_inametyrphostin-AG-1478  0.115920192   0.12512476  0.92643686
      differential_expression <- as.data.frame(coef(summary(LMM_output)))
      # crea la colonna drug dai nomi delle righe dei coefficienti dell'analisi "lmer"
      differential_expression$drug <- rownames(differential_expression)
      # elimina la prima riga contenente l'intercetta
      differential_expression <- differential_expression[2:nrow(differential_expression),]
      # elimina il prefisso pert_iname dal nome del farmaco nella colonna drug
      differential_expression$drug <- substr(differential_expression$drug, 11, nchar(differential_expression$drug))
      if (obj_is_na(gene_id)) {
        # lascia solo le colonne contenenti i dati necessari
        differential_expression <- differential_expression[, c("drug", "Estimate", "Std. Error", "t value")]
        # rinomina le colonne
        colnames(differential_expression) <- c("drug", "DE_log2_FC", "std.error", "t.value")
      }else {
        differential_expression$gene_id <- gene_id
        # lascia solo le colonne contenenti i dati necessari
        differential_expression <- differential_expression[, c("drug", "gene_id", "Estimate", "Std. Error", "t value")]
        # rinomina le colonne
        colnames(differential_expression) <- c("drug", "gene_id", "DE_log2_FC", "std.error", "t.value")
      }

      # elimina i nomi delle righe
      rownames(differential_expression) <- NULL
      # risultato finale:
      # drug                gene  Estimate      Std. Error  t value
      # tyrphostin-AG-1296  PAX8  -0.032810809  0.12512476  -0.26222474
      # tyrphostin-AG-1478  PAX8  0.115920192   0.12512476  0.92643686
      return(differential_expression)
    }
  )
)