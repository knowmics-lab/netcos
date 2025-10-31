
DGEMetanalysisMultiple <- R6Class(
  "DGEMetanalysisMultiple",
  public = list(
    initialize = function() {
      private$dgeMetanalysis <- DGEMetanalysis$new()
    },
    compute = function(drug_by_gene_DGE_A, drug_by_gene_DGE_B) {
      startTime <- Sys.time()
      dgrpLogger$log("start meta analysis")
      ### Combine 6h and 24h using meta-analysis: ###

      # calcola l'intersezione dei nomi dei farmaci tra gli insieme drugs_6h e drugs_24h.
      # merge fa la natural join
      drug_by_gene_DGE_A_B <- merge(x = drug_by_gene_DGE_A, y = drug_by_gene_DGE_B, by = c("gene_id", "drug"))
      # assegna i nomi alle colonne
      colnames(drug_by_gene_DGE_A_B)[3:12] <- c("DE_log2_FC_A", "std.error_A", "t.value_A", "p.value_A", "adj.p.value_A")
      colnames(drug_by_gene_DGE_A_B)[8:12] <- c("DE_log2_FC_B", "std.error_B", "t.value_B", "p.value_B", "adj.p.value_B")

      # in questo ciclo for viene fatta la meta analisi per calcolare l'effetto combinato di drugs_6h e drugs_24h
      total_drug_by_gene_dges <- nrow(drug_by_gene_DGE_A_B)
      metalalysis_A_B <- data.frame(matrix(NA, nrow = total_drug_by_gene_dges, ncol = 17))
      dgrpLogger$log(sprintf("dge common drugs: %s", total_drug_by_gene_dges))
      for (i in 1:total_drug_by_gene_dges) {
        metalalysis_A_B[i,] <- private$dgeMetanalysis$compute(drug_by_gene_DGE_A_B[i,])
      }
      colnames(metalalysis_A_B) <- c(colnames(drug_by_gene_DGE_A_B), c("DE_log2_FC_A_B", "std.error_A_B", "t.value_A_B", "p.value_A_B","adj.p.value_A_B"))
      # drug_by_gene_DGE_6h_24h_ma
      #    drug      gene   DE_log2_FC_6h    std.error_6h   t.value_6h      DE_log2_FC_24h    std.error_24h        t.value_24h      DE_log2_FC_6h_24h  std.error_6h_24h    t.value_6h_24h       p.value_6h_24h
      # 1  ABT-737   AARS   0.4461631309544  0.14302313904  0.155102125097  0.15510212509716  0.15479839964925124  1.0019620709813  0.3066570115074    0.1454057567166199  0.21089743517179400  0.0349467955686023
      # 2  ABT-751   AARS   0.0103021782160  0.05283095800  0.195002676573  -0.0534011868235  0.09920129528376088  -0.538311386668  -0.003773423983    0.0466304544571578  -0.0809218787845287  0.9355040802701673
      totalTime <- Sys.time() - startTime
      dgrpLogger$log(sprintf("end meta analysis, computation time: %s %s", totalTime, attr(totalTime, "units")))
      return(metalalysis_A_B)
    }
  ),
  private = list(
    dgeMetanalysis = NA
  )
)