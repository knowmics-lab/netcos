library(R6)

DiseaseSignatureMapper <- R6Class(
  "DiseaseSignatureMapper",
  public = list(
    map = function(disease_signature) {
      colnames(disease_signature)[1] <- "DE_log2_FC"
      disease_signature$gene <- rownames(disease_signature)
      disease_signature$DE_log2_FC_SE <- disease_signature$DE_log2_FC / disease_signature$t
      disease_signature <- disease_signature[, c("gene", "DE_log2_FC", "DE_log2_FC_SE", "t", "P.Value", "adj.P.Val")]
      colnames(disease_signature) <- c("gene", "DE_log2_FC", "std.error", "t.value", "p.value", "adj.p.value")
      disease_signature <- disease_signature[order(disease_signature$gene),]
      rownames(disease_signature) <- NULL
      return(disease_signature)
    }
  )
)