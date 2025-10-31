
VoomDGEMapper <- R6Class(
  "VoomDGEMapper",
  public = list(
    map = function(differential_expression) {
      differential_expression$gene_id <- rownames(differential_expression)
      differential_expression <- data.frame(differential_expression)
      differential_expression <- differential_expression[, c("gene_id", "logFC", "std.error", "t", "P.Value", "adj.P.Val")]
      colnames(differential_expression) <- c("gene_id", "DE_log2_FC", "std.error", "t.value", "p.value", "adj.p.value")
      differential_expression <- differential_expression[order(differential_expression$gene_id),]
      rownames(differential_expression) <- NULL
      return(differential_expression)
    }
  )
)
