# ----------------------------------------------------------------------------------------------
# Log2 transformation of GSE65000 FPKM matrix, using external_gene_name as final identifier
# ----------------------------------------------------------------------------------------------

# inputs
input_file  <- "D:/GroupA_University2026_Project/BulkRNA/GSE65000_hsa_fpkm_matrix.txt"
output_file <- "D:/GroupA_University2026_Project/BulkRNA/Old_GSE65000_hsa_fpkm_matrix_log2.txt"

id_column_to_keep_values_from <- "external_gene_id"   # the column whose values become the new gene IDs
id_column_name_in_output      <- "gene_symbol"          # name of the gene ID column in output file
columns_to_drop_always        <- c("ensembl_gene_id", "external_gene_id")  # these will be removed after use

data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Extract gene symbols
gene_symbols <- data[[id_column_to_keep_values_from]]

# Remove all non-numeric identifier columns no longer needed
cols_to_drop <- intersect(columns_to_drop_always, colnames(data))
all_non_numeric <- names(data)[!sapply(data, is.numeric)]
additional_drops <- setdiff(all_non_numeric, id_column_to_keep_values_from)
cols_to_drop <- unique(c(cols_to_drop, additional_drops))

if (length(cols_to_drop) > 0) {
  data <- data[, !(colnames(data) %in% cols_to_drop), drop = FALSE]
  cat("Dropped non-numeric columns:", paste(cols_to_drop, collapse = ", "), "\n")
}

expr_matrix <- data[, numeric_cols, drop = FALSE]

# Apply log2 transformation
log2_expr <- log2(expr_matrix)
log2_expr[log2_expr == -Inf] <- 0
log2_expr[is.na(log2_expr)] <- 0

output_df <- cbind(
  setNames(data.frame(gene_symbols, stringsAsFactors = FALSE), id_column_name_in_output),
  log2_expr
)

# Write the matrix
cat("Writing log2‑transformed data to", output_file, "...\n")
write.table(output_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)