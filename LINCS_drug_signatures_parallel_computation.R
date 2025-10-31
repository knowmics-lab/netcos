library(parallel)
source("modules/config.R")
source("modules/drug_signature/DrugSignatureInputArguments.R")
source("modules/drug_signature/splitter/GeneSymbolsVectorSplitter.R")
source("modules/drug_signature/builder/DrugSignatureInputChunksBuilder.R")
source("modules/TSRLogger.R")

# setup
drugSignatureInputArguments <- DrugSignatureInputArguments$new(drug_signature_cfg$gene_symbols_filename)
geneSymbolsVectorSplitter <- GeneSymbolsVectorSplitter$new()
drugSignatureInputChunksBuilder <- DrugSignatureInputChunksBuilder$new()
tsrLogger <- TSRLogger$new()

# computation
inputArguments <- drugSignatureInputArguments$get()
tsrLogger$log(sprintf("gene symbols file: %s", inputArguments$genes_filename))
perturbation_times_label <- paste0(paste(inputArguments$perturbation_times, collapse = "h, "), "h")
tsrLogger$log(sprintf("perturbation times: %s", perturbation_times_label))
BLAS_num_threads_label <- paste(inputArguments$BLAS_num_threads, collapse = ", ")
tsrLogger$log(sprintf("BLAS_num_threads: %s", BLAS_num_threads_label))
tsrLogger$log(sprintf("num_chunks: %s", inputArguments$num_chunks))
genes_symbols_chunks <- geneSymbolsVectorSplitter$split(inputArguments$genes_filename, inputArguments$num_chunks)
BLAS_num_threads <- inputArguments$BLAS_num_threads
num_process <- inputArguments$num_chunks * length(inputArguments$perturbation_times)
tsrLogger$log(sprintf("num processes (num_chunks * length(perturbation_times)): %s", num_process))

if (length(BLAS_num_threads) > 1 & num_process != length(BLAS_num_threads)) {
  stop("invalid input: bad BLAS_num_threads/num_chunks")
}

if (length(BLAS_num_threads) == 1 & inputArguments$num_chunks > 1) {
  BLAS_num_threads <- rep(BLAS_num_threads, num_process)
}

drug_signature_input_chunks <- drugSignatureInputChunksBuilder$build(genes_symbols_chunks, inputArguments$perturbation_times, BLAS_num_threads)

process_chunk <- function(chunk) {
  source("modules/config.R")
  delay <- (chunk$number - 1) * parallel_computation$delay_start_clusters
  print(sprintf("deleay before start: %s", delay))
  Sys.sleep(delay)
  source("modules/drug_signature/OverallDrugSignatureParallelMain.R")
  # setup
  overallDrugSignatureParallelMain <- OverallDrugSignatureParallelMain$new()

  # computation
  overallDrugSignatureParallelMain$compute(chunk$gene_symbols, chunk$perturbation_time, chunk$BLAS_num_threads)
}

cl <- makePSOCKcluster(num_process, outfile = "")
clusterApply(cl, drug_signature_input_chunks, process_chunk)
stopCluster(cl)
