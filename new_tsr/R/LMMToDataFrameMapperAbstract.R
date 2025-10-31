
LMMToDataFrameMapperAbstract <- R6Class(
  "LMMToDataFrameMapperAbstract",
  public = list(
    map = function(LMM_output, gene_symbol = NA) {
      stop("I'm an abstract method, implement me")
    }
  )
)