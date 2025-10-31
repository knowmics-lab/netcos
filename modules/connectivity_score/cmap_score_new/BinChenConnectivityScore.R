library(R6)
source("modules/connectivity_score/cmap_score_new/cmap_score_new.R")

BinChenConnectivityScore <- R6Class(
  "BinChenConnectivityScore",
  public = list(
    compute = function(down_regulated_genes, up_regulated_genes, drug_signature, alt_cs_distribution) {
      drug_signature$ids <- rownames(drug_signature)
      drug_signature$rank <- rank(-1 * drug_signature$estimate, ties.method = "random")

      cs <-
        cmap_score_new(
          up_regulated_genes,
          down_regulated_genes,
          drug_signature[, c("ids", "rank")]
        )

      output <- c(
        cs,
        sum(alt_cs_distribution < cs) / length(alt_cs_distribution)
      )

      names(output) <- c("connectivity_score", "connectivity_score_p.value")
      return(output)
    }

  )
)






