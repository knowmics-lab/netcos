
LINCSJuliaDrugDGEFactory <- R6Class(
  "LINCSJuliaDrugDGEFactory",
  public = list(
    initialize = function() {
      private$juliaSetuper <- JuliaSetuper$new()
    },
    create = function(BLAS_num_threads = NA) {
      private$juliaSetuper$setup(BLAS_num_threads)
      juliaLMM <- JuliaLMM$new(config$LINCSLMMFormula)
      juliaLMMToDataFrameMapper <- LINCSJuliaLMMToDataFrameMapper$new()
      drugDGE <- LMMDGEByGene$new(juliaLMM, juliaLMMToDataFrameMapper)
      return(drugDGE)
    }
  ),
  private = list(
    juliaSetuper = NA
  )
)