library(R6)
source("modules/drug_signature/julia/setuper/JuliaSetuper.R")
source("modules/drug_signature/julia/JuliaLMM.R")
source("modules/drug_signature/julia/mapper/JuliaLMMToDataFrameMapper.R")
source("modules/drug_signature/DrugSignature.R")

JuliaDrugSignatureFactory <- R6Class(
  "JuliaDrugSignatureFactory",
  public = list(
    initialize = function(experimentDataLoader, signature_dir = NA) {
      private$juliaSetuper <- JuliaSetuper$new()
      private$experimentDataLoader <- experimentDataLoader
      private$signature_dir <- signature_dir
    },
    create = function(BLAS_num_threads = NA) {
      private$juliaSetuper$setup(BLAS_num_threads)
      juliaLMM <- JuliaLMM$new()
      juliaLMMToDataFrameMapper <- JuliaLMMToDataFrameMapper$new()
      drugSignature <- DrugSignature$new(private$experimentDataLoader, juliaLMM, juliaLMMToDataFrameMapper, private$signature_dir)
      return(drugSignature)
    }
  ),
  private = list(
    juliaSetuper = NA,
    experimentDataLoader = NA,
    signature_dir = NA
  )
)