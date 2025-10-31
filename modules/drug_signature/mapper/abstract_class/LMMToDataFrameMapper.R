library(R6)

LMMToDataFrameMapper <- R6Class(
  "LMMToDataFrameMapper",
  public = list(
    map = function(LMM_output, gene) {
      stop("I'm an abstract method, implement me")
    }
  )
)