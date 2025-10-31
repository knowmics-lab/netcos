library(R6)
library(JuliaCall)

source("modules/drug_signature/mapper/abstract_class/LMMToDataFrameMapper.R")

JuliaLMMToDataFrameMapper <- R6Class(
  "JuliaLMMToDataFrameMapper",
  inherit = LMMToDataFrameMapper,
  public = list(
    map = function(LMM_output, gene) {
      LMM_table <- julia_call("coeftable", LMM_output)
      LMM_table <- JuliaCall::field(LMM_table, "cols")
      DE_log2_FC <- LMM_table[[1]][-1]
      std.error <- LMM_table[[2]][-1]
      t.value <- DE_log2_FC / std.error
      p.value <- LMM_table[[4]][-1]
      adj.p.value <- p.adjust(p.value, method = "BH")
      drugs <- julia_call("coefnames", LMM_output)[-1]
      drugs <- substr(drugs, 13, nchar(drugs))
      return(data.frame(
        drug = drugs,
        gene = gene,
        DE_log2_FC = DE_log2_FC,
        std.error = std.error,
        t.value = t.value,
        p.value = p.value,
        adj.p.value = adj.p.value
      ))
    }
  )
)