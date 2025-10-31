
LINCSJuliaLMMToDataFrameMapper <- R6Class(
  "LINCSJuliaLMMToDataFrameMapper",
  inherit = LMMToDataFrameMapperAbstract,
  public = list(
    map = function(LMM_output, gene_id = NA) {
      LMM_table <- julia_call("coeftable", LMM_output)
      LMM_table <- JuliaCall::field(LMM_table, "cols")
      DE_log2_FC <- LMM_table[[1]][-1]
      std.error <- LMM_table[[2]][-1]
      t.value <- DE_log2_FC / std.error
      p.value <- LMM_table[[4]][-1]
      adj.p.value <- p.adjust(p.value, method = "BH")
      drugs <- julia_call("coefnames", LMM_output)[-1]
      drugs <- substr(drugs, 13, nchar(drugs))
      if (obj_is_na(gene_id)) {
        dge <- data.frame(
          drug = drugs,
          DE_log2_FC = DE_log2_FC,
          std.error = std.error,
          t.value = t.value,
          p.value = p.value,
          adj.p.value = adj.p.value
        )
      }else {
        dge <- data.frame(
          drug = drugs,
          gene_id = gene_id,
          DE_log2_FC = DE_log2_FC,
          std.error = std.error,
          t.value = t.value,
          p.value = p.value,
          adj.p.value = adj.p.value
        )
      }
      return(dge)
    }
  )
)
