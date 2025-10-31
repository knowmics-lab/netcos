IlluminaDichotomicVoomDrugRepurpose <- R6Class(
  "IlluminaDichotomicVoomDrugRepurpose",
  public = list(
    compute = function(rna_seq_filename, sample_01_map, disease_name,
                       n_most_significant_genes, drug_dge_dir,
                       drugs, drugs_genes, drug_gde_t_value_column_name = "t.value",
                       random_connectivity_score_permutations = 10^5,
                       drug_perturbation_time = NA, parallel_computation = F,
                       filter_by_protein_coding = F
    ) {
      illuminaDichotomicVoomDGE <- IlluminaDichotomicVoomDGE$new()
      drugSignatureLoaderByDrugName <- DrugSignatureLoaderByDrugName$new(drug_dge_dir, drug_gde_t_value_column_name)
      overallDiseaseDrugConnectivityScore <- OverallDiseaseDrugConnectivityScore$new(drugSignatureLoaderByDrugName)
      startTime <- Sys.time()
      dgrpLogger$log("start drug repurposing computation")
      disease_dge <- illuminaDichotomicVoomDGE$compute(rna_seq_filename, sample_01_map, disease_name, filter_by_protein_coding)
      connectivity_score <- overallDiseaseDrugConnectivityScore$compute(
        disease_dge, drugs, drugs_genes, n_most_significant_genes, random_connectivity_score_permutations, disease_name, drug_perturbation_time, parallel_computation)
      totalTime <- Sys.time() - startTime
      dgrpLogger$log(sprintf("end drug repurposing computation, time: %s %s", totalTime, attr(totalTime, "units")))
      return(connectivity_score)
    }
  )
)