# -------------------------------------------------------------------
# Read matrix and list samples not expressing ACTB or GAPDH
# -------------------------------------------------------------------

# Inputs
matrix_file <- "2026_UoL_prj_gene_matrix_f.tsv"
threshold   <- -1               # log2(0.5) – samples with value ≤ this (or NA) are "not expressed"
hk_genes    <- c("ACTB", "GAPDH")

log2_mat <- read.delim(matrix_file, row.names = 1, check.names = FALSE)

# Ensure the housekeeping genes are present
missing_in_matrix <- setdiff(hk_genes, rownames(log2_mat))
if (length(missing_in_matrix) > 0) {
  stop("Gene(s) not found in the matrix: ", paste(missing_in_matrix, collapse = ", "))
}

actb_vec  <- as.numeric(log2_mat["ACTB", ])
gapdh_vec <- as.numeric(log2_mat["GAPDH", ])

# Identify samples where expression is missing or ≤ threshold
actb_fail  <- is.na(actb_vec)  | actb_vec  <= threshold
gapdh_fail <- is.na(gapdh_vec) | gapdh_vec <= threshold

failed_samples <- colnames(log2_mat)[actb_fail | gapdh_fail]

# Print results
if (length(failed_samples) == 0) {
  cat("All samples express both housekeeping genes (log2(FPKM) > -1).\n")
} else {
  cat("Samples not expressing ACTB and/or GAPDH:\n")
  for (samp in failed_samples) {
    reason <- c()
    if (actb_fail[samp])  reason <- c(reason, "ACTB")
    if (gapdh_fail[samp]) reason <- c(reason, "GAPDH")
    cat(sprintf("  %s – %s\n", samp, paste(reason, collapse = ", ")))
  }
}