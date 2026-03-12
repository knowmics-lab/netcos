#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(metafor)
})

stopf <- function(...) stop(sprintf(...), call. = FALSE)
warnf <- function(...) warning(sprintf(...), call. = FALSE)

standardize_signature_df <- function(x, time_label, drug_name = NA_character_) {
  if (!is.data.frame(x)) {
    stopf("[%s] RDS object is not a data.frame (got: %s)", time_label, class(x)[1])
  }
  
  required <- c("drug", "gene_id", "DE_log2_FC", "std.error", "t.value", "p.value", "adj.p.value")
  missing <- setdiff(required, colnames(x))
  if (length(missing) > 0) {
    stopf("[%s] Missing required columns: %s", time_label, paste(missing, collapse = ", "))
  }
  
  out <- data.frame(
    drug       = as.character(x$drug),
    gene_id    = as.character(x$gene_id),
    DE_log2_FC = as.numeric(x$DE_log2_FC),
    std.error  = as.numeric(x$std.error),
    t.value    = as.numeric(x$t.value),
    p.value    = as.numeric(x$p.value),
    adj.p.value = as.numeric(x$adj.p.value),
    stringsAsFactors = FALSE
  )
  
  # Drop missing/empty gene_id
  out <- out[!is.na(out$gene_id) & nzchar(out$gene_id), , drop = FALSE]
  out
}

compute_metanalysis_by_gene <- function(df6, df24, ma_obj) {
  # Merge by gene_id and drug
  m <- merge(df6, df24, by = c("gene_id", "drug"))
  if (nrow(m) == 0) {
    return(m)  # returns empty data.frame with correct columns
  }
  
  # Rename columns exactly like TSR expects
  colnames(m) <- c(
    "gene_id","drug",
    "DE_log2_FC_6h","std.error_6h","t.value_6h","p.value_6h","adj.p.value_6h",
    "DE_log2_FC_24h","std.error_24h","t.value_24h","p.value_24h","adj.p.value_24h"
  )
  
  total <- nrow(m)
  out <- data.frame(matrix(NA, nrow = total, ncol = ncol(m) + 4))
  
  for (i in seq_len(total)) {
    out[i, ] <- ma_obj$compute(m[i, , drop = FALSE])
  }
  
  colnames(out) <- c(
    colnames(m),
    "DE_log2_FC_6h_24h","std.error_6h_24h","t.value_6h_24h","p.value_6h_24h"
  )
  
  out
}


# -----------------------------
# Main
# -----------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat(
    "Usage:\n",
    "  Rscript metanalysis_6h_24h_batch.R <base_dir> <cell_line_name> [<netcos_root>]\n\n",
    "Expected folder structure inside <base_dir>:\n",
    "  cell_line_6h/\n",
    "  cell_line_24h/\n",
    "Output will be created as:\n",
    "  cell_line_metanalysis/\n\n",
    "Example:\n",
    "  Rscript metanalysis_6h_24h_batch.R data/A549 A549 .\n",
    sep = ""
  )
  quit(status = 1)
}

base_dir <- normalizePath(args[[1]], mustWork = TRUE)
cell_line_name <- args[[2]]
netcos_root <- if (length(args) >= 3) args[[3]] else "."

# --- NEW: robust folder resolver ---
resolve_time_dirs <- function(base_dir, cell_line_name) {
  candidates_6h <- c(
    file.path(base_dir, paste0(cell_line_name, "_6h")),
    file.path(base_dir, paste0(cell_line_name, "-6h")),
    file.path(base_dir, paste0(cell_line_name, "6h")),
    file.path(base_dir, cell_line_name, "6h"),
    file.path(base_dir, cell_line_name, "cell_line_6h"),
    file.path(base_dir, "cell_line_6h")
  )
  candidates_24h <- c(
    file.path(base_dir, paste0(cell_line_name, "_24h")),
    file.path(base_dir, paste0(cell_line_name, "-24h")),
    file.path(base_dir, paste0(cell_line_name, "24h")),
    file.path(base_dir, cell_line_name, "24h"),
    file.path(base_dir, cell_line_name, "cell_line_24h"),
    file.path(base_dir, "cell_line_24h")
  )
  
  pick_existing <- function(paths) {
    paths <- normalizePath(paths, mustWork = FALSE)
    paths[dir.exists(paths)][1]
  }
  
  dir_6h <- pick_existing(candidates_6h)
  dir_24h <- pick_existing(candidates_24h)
  
  # fallback: search under base_dir for dirs containing both cell line and time token
  if (is.na(dir_6h) || is.na(dir_24h)) {
    all_dirs <- unique(list.dirs(base_dir, recursive = TRUE, full.names = TRUE))
    
    match_time_dir <- function(token) {
      hits <- all_dirs[
        grepl(cell_line_name, basename(all_dirs), ignore.case = TRUE) &
          grepl(token, basename(all_dirs), ignore.case = TRUE)
      ]
      if (length(hits) == 0) return(NA_character_)
      hits[[1]]
    }
    
    if (is.na(dir_6h))  dir_6h  <- match_time_dir("6h")
    if (is.na(dir_24h)) dir_24h <- match_time_dir("24h")
  }
  
  if (is.na(dir_6h) || is.na(dir_24h)) {
    stop(sprintf(
      "Could not auto-detect 6h/24h folders under: %s\nFound:\n  6h = %s\n  24h = %s\n",
      base_dir,
      ifelse(is.na(dir_6h), "NOT FOUND", dir_6h),
      ifelse(is.na(dir_24h), "NOT FOUND", dir_24h)
    ), call. = FALSE)
  }
  
  list(dir_6h = dir_6h, dir_24h = dir_24h)
}

dirs <- resolve_time_dirs(base_dir, cell_line_name)
dir_6h <- dirs$dir_6h
dir_24h <- dirs$dir_24h

# output folder (created automatically)
out_dir <- file.path(base_dir, paste0(cell_line_name, "_metanalysis"))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message(sprintf("[%s] Using 6h folder:  %s", cell_line_name, dir_6h))
message(sprintf("[%s] Using 24h folder: %s", cell_line_name, dir_24h))
message(sprintf("[%s] Output folder:    %s", cell_line_name, out_dir))



# Source your pipeline class
metanalysis_r_path <- file.path(netcos_root, "netcos", "modules", "drug_signature", "Metanalysis6h24hDrugSignature.R")
if (!file.exists(metanalysis_r_path)) {
  metanalysis_r_path2 <- file.path(netcos_root, "modules", "drug_signature", "Metanalysis6h24hDrugSignature.R")
  if (!file.exists(metanalysis_r_path2)) {
    stopf("Cannot find Metanalysis6h24hDrugSignature.R at:\n- %s\n- %s", metanalysis_r_path, metanalysis_r_path2)
  }
  metanalysis_r_path <- metanalysis_r_path2
}
source(metanalysis_r_path)

ma <- Metanalysis6h24hDrugSignature$new()

files_6h  <- list.files(dir_6h,  pattern = "\\.Rds$", full.names = TRUE)
files_24h <- list.files(dir_24h, pattern = "\\.Rds$", full.names = TRUE)

b6  <- basename(files_6h)
b24 <- basename(files_24h)

common <- intersect(b6, b24)
if (length(common) == 0) stopf("[%s] No common .Rds filenames found between:\n- %s\n- %s", cell_line_name, dir_6h, dir_24h)

only_6h  <- setdiff(b6, b24)
only_24h <- setdiff(b24, b6)
if (length(only_6h) > 0)  warnf("[%s] Files only in 6h folder (skipped): %s",  cell_line_name, paste(head(only_6h, 20), collapse=", "))
if (length(only_24h) > 0) warnf("[%s] Files only in 24h folder (skipped): %s", cell_line_name, paste(head(only_24h, 20), collapse=", "))

for (bn in common) {
  drug_name <- sub("\\.Rds$", "", bn)
  
  f6  <- files_6h[match(bn, b6)]
  f24 <- files_24h[match(bn, b24)]
  
  message(sprintf("[%s] Processing %s", cell_line_name, drug_name))
  
  x6  <- readRDS(f6)
  x24 <- readRDS(f24)
  
  s6  <- standardize_signature_df(x6,  "6h",  drug_name = drug_name)
  s24 <- standardize_signature_df(x24, "24h", drug_name = drug_name)
  
  # merged <- merge(
  #   s6, s24,
  #   by = "gene_id",
  #   suffixes = c("_6h", "_24h"),
  #   all = FALSE
  # )
  # 
  # 
  # if (nrow(merged) == 0) {
  #   warnf("[%s] %s: No overlapping genes between 6h and 24h; skipping.", cell_line_name, drug_name)
  #   next
  # }
  # 
  # if ("drug_6h" %in% colnames(merged)) {
  #   merged$drug <- merged$drug_6h
  # } else if ("drug_24h" %in% colnames(merged)) {
  #   merged$drug <- merged$drug_24h
  # } else {
  #   merged$drug <- drug_name
  # }
  
  # OPTIONAL: keep only rows where both timepoints have finite FC and positive SE
  s6_ok  <- is.finite(s6$DE_log2_FC) & is.finite(s6$std.error) & s6$std.error > 0
  s24_ok <- is.finite(s24$DE_log2_FC) & is.finite(s24$std.error) & s24$std.error > 0
  s6_f  <- s6[s6_ok, , drop = FALSE]
  s24_f <- s24[s24_ok, , drop = FALSE]
  
  if (nrow(s6_f) == 0 || nrow(s24_f) == 0) {
    warnf("[%s] %s: No valid rows after filtering; skipping.", cell_line_name, drug_name)
    next
  }
  
  result <- compute_metanalysis_by_gene(s6_f, s24_f, ma)
  
  saveRDS(result, file.path(out_dir, bn))
  
  # out <- merged
  # out$DE_log2_FC_6h_24h <- NA_real_
  # out$std.error_6h_24h  <- NA_real_
  # out$t.value_6h_24h    <- NA_real_
  # out$p.value_6h_24h    <- NA_real_
  # 
  # valid <- is.finite(out$DE_log2_FC_6h) & is.finite(out$DE_log2_FC_24h) &
  #   is.finite(out$std.error_6h)  & is.finite(out$std.error_24h)  &
  #   out$std.error_6h  > 0 & out$std.error_24h > 0
  # 
  # if (any(valid)) {
  #   #tmp <- out[valid, , drop = FALSE]
  #   #tmp <- ma$compute(tmp)
  #   
  #   out[valid, c("DE_log2_FC_6h_24h", "std.error_6h_24h", "t.value_6h_24h", "p.value_6h_24h")] <-
  #     tmp[, c("DE_log2_FC_6h_24h", "std.error_6h_24h", "t.value_6h_24h", "p.value_6h_24h")]
  # } else {
  #   warnf("[%s] %s: No valid rows for rma() (check std.error columns). Writing NA-only metanalysis columns.",
  #         cell_line_name, drug_name)
  # }
  # 
  # saveRDS(out, file.path(out_dir, bn))
}

message(sprintf("[%s] Done. Metanalysis saved in: %s", cell_line_name, out_dir))
