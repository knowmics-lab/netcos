
DGEMetanalysis <- R6Class(
  "DGEMetanalysis",
  public = list(
    compute = function(dge) {
      ma_output <- rma(yi = c(dge$DE_log2_FC_A, dge$DE_log2_FC_B), sei = c(dge$std.error_A, dge$std.error_B))
      dge[c("DE_log2_FC_A_B", "std.error_A_B", "t.value_A_B", "p.value_A_B")] <-
        c(ma_output$b[, 1], ma_output$se, ma_output$zval, ma_output$pval)
      dge["adj.p.value_A_B"] <- p.adjust(dge["p.value_A_B"], method = "BH")
      return(dge)
    }
  )
)