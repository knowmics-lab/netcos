
ConnectivityScorePValue <- R6Class(
  "ConnectivityScorePValue",
  public = list(
    compute = function(connectivity_score, random_cs_distribution) {
      return(sum(random_cs_distribution < connectivity_score) / length(random_cs_distribution))
    }
  )
)






