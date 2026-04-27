# -------------------------------
# Filter protein-coding genes
# -------------------------------

# Inputs
tsv_data <- read.delim(
  "D:/GroupA_University2026_Project/PCA_Output/Filtered_Expression_matrix.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

genes_raw <- read.delim(
  "D:/GroupA_University2026_Project/PCA_Output/2024_Coding_genes_list.txt",
  header = FALSE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Convert entire file into a simple vector
genes_to_keep <- as.vector(as.matrix(genes_raw))

# Clean gene names
genes_to_keep <- trimws(genes_to_keep)
genes_to_keep <- gsub('"', '', genes_to_keep)
genes_to_keep <- genes_to_keep[!is.na(genes_to_keep)]
genes_to_keep <- genes_to_keep[genes_to_keep != ""]

tsv_data[,1] <- trimws(tsv_data[,1])
tsv_data[,1] <- gsub('"', '', tsv_data[,1])

# Filter genes
filtered_data <- tsv_data[tsv_data[,1] %in% genes_to_keep, ]

# Diagnostics
cat("Total genes in matrix:", nrow(tsv_data), "\n")
cat("Genes in list:", length(genes_to_keep), "\n")
cat("Matched genes:", nrow(filtered_data), "\n")

write.table(
  filtered_data,
  "D:/GroupA_University2026_Project/PCA_Output/Filtered_Expression_matrix.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)