
LMMDGEByGene <- R6Class(
  "LMMDGEByGene",
  public = list(
    initialize = function(lmm, lmmToDataFrameMapper) {
      if (!"LMMAbstract" %in% class(lmm))
        stop("the lmm instance must by of type LMMAbstract")
      if (!"LMMToDataFrameMapperAbstract" %in% class(lmmToDataFrameMapper))
        stop("the lmmToDataFrameMapper instance must by of type LMMToDataFrameMapperAbstract")
      private$lmm <- lmm
      private$lmmToDataFrameMapper <- lmmToDataFrameMapper
    },
    compute = function(rna_data, gene_symbol = NA) {
      dgrpLogger$log(sprintf("start differential gene expression computation %s", gene_symbol))
      startTime <- Sys.time()
      # calcola il modello lineare misto:
      # vadiabile dipendente: "expr", espressione genica
      # effetto fisso: "pert_iname", composto chimico (farmaco)
      # effetti random: cell_id, rna_plate
      # LMM_output <- lmer(gene_expressions ~ pert_iname + (1 | cell_id) + (1 | rna_plate), data = experiments_metadata, control = lmerControl(calc.derivs = FALSE))

      differentialExpression <- private$lmm$compute(rna_data)

      differentialExpression <- private$lmmToDataFrameMapper$map(differentialExpression, gene_symbol)
      # risultato finale:
      # drug                gene  Estimate      Std. Error  t value
      # tyrphostin-AG-1296  PAX8  -0.032810809  0.12512476  -0.26222474
      # tyrphostin-AG-1478  PAX8  0.115920192   0.12512476  0.92643686
      totalTime <- Sys.time() - startTime
      dgrpLogger$log(sprintf("end gene differential expression computation %s, time: %s %s", gene_symbol, totalTime, attr(totalTime, "units")))
      # force garbage collection to prevent out of memory
      gc()
      return(differentialExpression)
    }
  ),
  private = list(
    lmm = NA,
    lmmToDataFrameMapper = NA
  )
)
