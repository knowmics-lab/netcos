
LINCSGCTXDataRowLoader <- R6Class(
  "LINCSGCTXDataRowLoader",
  public = list(
    initialize = function(experiments_data_filename) {
      private$experiments_data_filename <- experiments_data_filename
    },
    load = function(gene_ids) {
      startTime <- Sys.time()
      gene_expressions <-
        parse_gctx(
          fname = private$experiments_data_filename,
          rid = gene_ids
        )@mat
      totalTime <- Sys.time() - startTime
      dgrpLogger$log(sprintf("end parsing gctx: %s %s", totalTime, attr(totalTime, "units")))
      return(gene_expressions)
    }
  ),
  private = list(
    experiments_data_filename = NA
  )
)
