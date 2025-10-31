
DGEMergeByGeneSync <- R6Class(
  "DGEMergeByGeneSync",
  public = list(
    initialize = function() {
      private$filenameBuilder <- FilenameBuilder$new()
    },
    merge = function(gene_list, group_list, input_dge_dir, output_group_dge_dir, group_field_name, dge_filename_pattern = "#id.Rds") {
      if (!private$filenameBuilder$valid_filename_pattern(dge_filename_pattern)) {
        stop("dge_filename_pattern should contain \"#id\" patten")
      }
      dgrpLogger$log("starting differential gene expressions reading")
      processorCores$initCores()
      total_genes <- length(gene_list)
      merged_dges <- data.frame()
      for (i in 1:total_genes) {
        filename <- private$filenameBuilder$build(dge_filename_pattern, gene_list[i])
        dgrpLogger$log(paste0(i, " - ", filename))
        filename <- paste0(input_dge_dir, filename)
        dge <- readRDS(filename)
        merged_dges <- rbind(merged_dges, subset(dge, dge[[group_field_name]] %in% group_list))
      }
      dgrpLogger$log("ending differential gene expressions reading")
      dgrpLogger$log("starting differential gene expressions group saving")
      total_groups <- length(group_list)
      for (i in 1:total_groups) {
        selected_group <- subset(merged_dges, merged_dges[[group_field_name]] %in% group_list[i])
        filename <- paste0(output_group_dge_dir, group_list[i], ".Rds")
        rownames(selected_group) <- NULL
        saveRDS(selected_group, filename)
        dgrpLogger$log(paste0(i, " - ", group_list[i]))
      }
      dgrpLogger$log("ending differential gene expressions group saving")
    }
  ),
  private = list(
    filenameBuilder = NA
  )
)
