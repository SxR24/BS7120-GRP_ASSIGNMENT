# -----------------------------------------------------------
# Combine Cufflinks files into a matrix and calculates log2(FPKM)
# -----------------------------------------------------------

dir.create("D:/GroupA_University2026_Project/PCA_Output", showWarnings = FALSE, recursive = TRUE)

# each subfolder named by sample ID and containing genes.fpkm_tracking
files <- list.files(
  path = "Cufflinks",
  pattern = "^genes\\.fpkm_tracking$",
  recursive = TRUE,
  full.names = TRUE
)

# List to hold a two-column data.frame for each sample (gene, expression)
expr_list <- list()

for (f in files) {
  sample_id <- basename(dirname(f))
  df <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  
  # gene identifier column (Use first one)
  gene_col <- intersect(c("gene_short_name", "gene_id", "tracking_id"), names(df))[1]
  
  # Convert FPKM to numeric. Log2(FPKM) (0 if <= 0 to avoid -inf)
  fpkm_num <- suppressWarnings(as.numeric(df$FPKM))
  log2_fpkm <- ifelse(!is.na(fpkm_num) & fpkm_num > 0, log2(fpkm_num), 0)
  
  tmp <- data.frame(
    gene = df[[gene_col]],
    value = log2_fpkm,
    stringsAsFactors = FALSE
  )
  tmp <- aggregate(value ~ gene, data = tmp, FUN = max)
  colnames(tmp)[2] <- sample_id
  expr_list[[sample_id]] <- tmp
}

# Merge all samples into a matrix
expr_mat_df <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), expr_list)
rownames(expr_mat_df) <- expr_mat_df$gene
expr_mat_df$gene <- NULL
expr_mat_df[is.na(expr_mat_df)] <- 0
expr_mat <- as.matrix(expr_mat_df)

write.table(
  cbind(gene = rownames(expr_mat), expr_mat),
  file = file.path("PCA_Output", "Expression_matrix_log2FPKM.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)