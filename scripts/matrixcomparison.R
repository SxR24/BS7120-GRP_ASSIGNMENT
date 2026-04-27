# =============================================================================
# Element‑wise comparison of original and papers expression matrices
# =============================================================================

library(readr)
library(dplyr)
library(tibble)


# Inputs
mat1_raw <- read_tsv("GSE75140_hOrg.fetal.master.data.frame.txt", 
                     show_col_types = FALSE) %>%
  as.data.frame()

mat2_raw <- read_tsv("2026_UoL_prj_gene_matrix_f.tsv", 
                     show_col_types = FALSE) %>%
  as.data.frame()

# Uniform orientation: rows = genes, columns = samples
rownames(mat1_raw) <- mat1_raw[, 1]
mat1_raw <- mat1_raw[, -1]
mat1 <- t(as.matrix(mat1_raw))   

rownames(mat2_raw) <- mat2_raw[, 1]
mat2_raw <- mat2_raw[, -1]
mat2 <- as.matrix(mat2_raw)



# Align by common genes and common samples
all_genes <- union(rownames(mat1), rownames(mat2))
all_samples <- union(colnames(mat1), colnames(mat2))

# Create empty matrices filled with NA
mat1_full <- matrix(NA_real_,
                    nrow = length(all_genes),
                    ncol = length(all_samples),
                    dimnames = list(all_genes, all_samples))
mat2_full <- mat1_full   

genes1 <- intersect(rownames(mat1), all_genes)
samples1 <- intersect(colnames(mat1), all_samples)
mat1_full[genes1, samples1] <- mat1[genes1, samples1]

genes2 <- intersect(rownames(mat2), all_genes)
samples2 <- intersect(colnames(mat2), all_samples)
mat2_full[genes2, samples2] <- mat2[genes2, samples2]

mat1_sub <- mat1_full
mat2_sub <- mat2_full

fc_threshold <- 2   

# Flatten each matrix into vectors 
vals1 <- as.vector(mat1_sub)
vals2 <- as.vector(mat2_sub)

keep <- is.finite(vals1) & is.finite(vals2) & vals1 > 0 & vals2 > 0
vals1 <- vals1[keep]
vals2 <- vals2[keep]

total_pairs <- length(vals1)
cat("Valid (sample, gene) pairs compared:", total_pairs, "\n")

if (total_pairs == 0) stop("No valid numeric pairs.")

# Calculate fold‑change: always ≥ 1 (max/min)
ratio <- pmax(vals1 / vals2, vals2 / vals1)

# Count pairs within the threshold
concordant <- sum(ratio <= fc_threshold)
percent_concordant <- (concordant / total_pairs) * 100

cat(sprintf("Pairs within %.1f‑fold: %d / %d (%.2f%%)\n",
            fc_threshold, concordant, total_pairs, percent_concordant))