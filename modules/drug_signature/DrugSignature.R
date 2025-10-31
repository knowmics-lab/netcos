library(R6)
source("modules/config.R")
source("modules/TSRLogger.R")
source("modules/obj_is_na.R")

DrugSignature <- R6Class(
  "DrugSignature",
  public = list(
    initialize = function(experimentDataLoader, lmm, lmmToDataFrameMapper, signature_dir = NA) {
      if (!"ExperimentDataLoader" %in% class(experimentDataLoader))
        stop("the experimentDataLoader instance must by of type ExperimentDataLoader")
      if (!"LMM" %in% class(lmm))
        stop("the lmm instance must by of type LMM")
      if (!"LMMToDataFrameMapper" %in% class(lmmToDataFrameMapper))
        stop("the lmmToDataFrameMapper instance must by of type LMMToDataFrameMapper")
      private$experimentDataLoader <- experimentDataLoader
      private$lmm <- lmm
      private$lmmToDataFrameMapper <- lmmToDataFrameMapper
      private$tsrLogger <- TSRLogger$new()
      if (obj_is_na(signature_dir)) {
        private$signature_dir <- paste0(drug_signature_cfg$signature_base_dir, drug_signature_cfg$genes_signature_dir)
      }else {
        private$signature_dir <- signature_dir
      }
    },
    compute = function(experiments_env, gene_id, gene_symbol, computation_number) {
      pert_time_hours <- experiments_env$experiments_meta_data$pert_time[1]
      filename <- paste0(private$signature_dir, gene_symbol, "_", gene_id, "_", pert_time_hours, "h.Rds")
      if (drug_signature_cfg$skip_already_computed_genes & file.exists(filename)) {
        private$tsrLogger$log(sprintf("skipping gene %s, purtubation time hours: %s, number: %s", gene_symbol, pert_time_hours, computation_number))
        return(NA)
      }
      private$tsrLogger$log(sprintf("start computation gene %s, purtubation time hours: %s, number: %s", gene_symbol, pert_time_hours, computation_number))
      gene_expressions <- private$experimentDataLoader$load(gene_id, experiments_env, computation_number)
      startTime <- Sys.time()

      private$tsrLogger$log(sprintf("start %s computation number: %s", private$lmm$algorithmName(), computation_number))
      # calcola il modello lineare misto:
      # vadiabile dipendente: "expr", espressione genica
      # effetto fisso: "pert_iname", composto chimico (farmaco)
      # effetti random: cell_id, rna_plate
      # LMM_output <- lmer(gene_expressions ~ pert_iname + (1 | cell_id) + (1 | rna_plate), data = experiments_env$experiments_meta_data, control = lmerControl(calc.derivs = FALSE))

      LMM_output <- private$lmm$compute(experiments_env, gene_expressions)
      #julia_command("fit!(LinearMixedModel(@formula(gene_expressions ~ pert_iname + (1 | cell_id) + (1 | rna_plate)), experiments_meta_data));")

      gene_DE <- private$lmmToDataFrameMapper$map(LMM_output, gene_symbol)
      # risultato finale:
      # drug                gene  Estimate      Std. Error  t value
      # tyrphostin-AG-1296  PAX8  -0.032810809  0.12512476  -0.26222474
      # tyrphostin-AG-1478  PAX8  0.115920192   0.12512476  0.92643686
      totalTime <- Sys.time() - startTime
      private$tsrLogger$log(sprintf("end computation gene %s, purtubation time hours: %s, number: %s, time: %s %s", gene_symbol, pert_time_hours, computation_number, totalTime, attr(totalTime, "units")))
      remove(LMM_output)
      # force garbage collection
      gc()
      #private$tsrLogger$log("garbage collector julia")
      #julia_command("GC.gc(true)")
      #julia_command("GC.gc(false)")
      saveRDS(gene_DE, file = filename)
      return(gene_DE)
    }
  ),
  private = list(
    experimentDataLoader = NA,
    lmm = NA,
    lmmToDataFrameMapper = NA,
    tsrLogger = NA,
    signature_dir = NA
  )
)
