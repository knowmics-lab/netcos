library(R6)

ConnectivityScoreBuilder <- R6Class(
  "ConnectivityScoreBuilder",
  public = list(
    build = function(disease_name, drugs, connectivity_score_matrix, gene_selection, perturbation_time) {
      drugs_connectivity_score <- data.frame(
        disease = disease_name,
        drug = drugs[connectivity_score_matrix[, 1]],
        gene_selection = gene_selection,
        perturbation_time = perturbation_time,
        connectivity_score = connectivity_score_matrix[, 2],
        p.value = connectivity_score_matrix[, 3]
      )
      rownames(drugs_connectivity_score) <- NULL
      return(drugs_connectivity_score)
      #     disease   drug      gene_selection        perturbation_time   connectivity_score   connectivity_score_p.value
      # 1   als       ABT-737   50_most_significant   6h                 0.013243241          0.0034293042000
      # 2   als       ABT-751   50_most_significant   6h                 -0.59877771          0.0048975924999
      # 3   als       afatinib  50_most_significant   6h                 -0.12349400          0.9004344124111
    }
  )
)