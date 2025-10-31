
RnaSeqSampleMapBuilder <- R6Class(
  "RnaSeqSampleMapBuilder",
  public = list(
    build = function(sample_01_map, test_sample_name) {
      # sample_selection_01_map = 001100X01 => 0 = test sample, 1 = control sample, X = excluded sample
      sample_vector_map <- strsplit(sample_01_map, split = "")[[1]]
      sample_positions <- which(sample_vector_map != "X")
      sample_vector_map <- subset(sample_vector_map, sample_vector_map != "X")
      sample_vector_map[which(sample_vector_map == '0')] <- test_sample_name
      sample_vector_map[which(sample_vector_map == '1')] <- "Control"
      samples <- factor(sample_vector_map, levels = c("Control", test_sample_name))
      return(
        list(
          samples = samples, # example: Control, Control, Test_sample, Comtol, ...
          sample_positions = sample_positions # 1, 3, 4, 8
        )
      )
    }
  ),
  private = list(
    testSampleControlGroupsFactor = NA
  )
)

