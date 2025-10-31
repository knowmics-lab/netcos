library(cmapR)

load_LINCS_experiment_data_rows <- function(experiments_env, experiments_data_filename, gene_ids) {
  startTime <- Sys.time()
  gene_expressions <-
    parse_gctx(
      fname = experiments_data_filename,
      rid = gene_ids,
      cid = experiments_env$experiments_meta_data$inst_id
    )@mat
  totalTime <- Sys.time() - startTime
  print(sprintf("end parsing gctx: %s %s", totalTime, attr(totalTime, "units")))
  return(gene_expressions)
}
