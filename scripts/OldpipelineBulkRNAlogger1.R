# ----------------------------------------------------------------------------------------------
# Log2 transformation of GSE38805 FPKM matrix, using external_gene_name as final identifier
# ----------------------------------------------------------------------------------------------

# Inputs
input_expr_file <- "D:/GroupA_University2026_Project/BulkRNA/GSE38805_human_FPKM.txt"
gtf_file        <- "D:/GroupA_University2026_Project/BulkRNA/Homo_sapiens.GRCh38.115.gtf"
output_file     <- "D:/GroupA_University2026_Project/BulkRNA/Old_GSE38805_human_FPKM_matrix_log2_symbols.txt"

expr <- read.table(input_expr_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Read GTF for gene_id to convert gene_name

gtf <- read.table(gtf_file, sep = "\t", header = FALSE, 
                  comment.char = "#", quote = "", stringsAsFactors = FALSE,
                  col.names = c("seqname", "source", "feature", "start", "end",
                                "score", "strand", "frame", "attributes"))
cat("  GTF loaded:", nrow(gtf), "total lines\n")

gtf_genes <- gtf[gtf$feature == "gene", ]
cat("  Gene entries in GTF:", nrow(gtf_genes), "\n")

# Extract gene_id and gene_name from attributes column
gene_ids_raw   <- sub('.*gene_id "([^"]+)".*', '\\1', gtf_genes$attributes)
gene_names     <- sub('.*gene_name "([^"]+)".*', '\\1', gtf_genes$attributes)

gene_ids_clean <- sub("\\.[0-9]+$", "", gene_ids_raw)

# Create mapping table
mapping <- data.frame(
  gene_id   = gene_ids_clean,
  gene_name = gene_names,
  stringsAsFactors = FALSE
)

mapping <- mapping[gene_ids_raw != gtf_genes$attributes & gene_names != gtf_genes$attributes, ]
mapping <- mapping[!duplicated(mapping$gene_id), ]
cat("  Unique gene_id entries in mapping:", nrow(mapping), "\n")

# Prepare expression matrix gene IDs
ensembl_ids_raw <- rownames(expr)
ensembl_ids_clean <- sub("\\.[0-9]+$", "", ensembl_ids_raw)

# Map to gene symbols
matched_idx <- match(ensembl_ids_clean, mapping$gene_id)
mapped_symbols <- mapping$gene_name[matched_idx]

n_mapped <- sum(!is.na(mapped_symbols))
n_missing <- sum(is.na(mapped_symbols))
cat("  Mapped", n_mapped, "genes to symbols (", n_missing, "missing)\n")

# For unmapped, keep original Ensembl Gene ID as identifier
mapped_symbols[is.na(mapped_symbols)] <- ensembl_ids_raw[is.na(mapped_symbols)]

# Handle duplicate gene symbols
dup_symbols <- duplicated(mapped_symbols) | duplicated(mapped_symbols, fromLast = TRUE)
if (any(dup_symbols)) {
  cat("  Found", sum(dup_symbols), "duplicate gene symbols. Making unique...\n")
  for (sym in unique(mapped_symbols[dup_symbols])) {
    idx <- which(mapped_symbols == sym)
    mapped_symbols[idx] <- paste0(sym, "_", ensembl_ids_raw[idx])
  }
}

# Replace row names with gene symbols
rownames(expr) <- mapped_symbols

# Apply log2 transformation
cat("Applying log2 transformation...\n")
log2_expr <- log2(expr)
log2_expr[log2_expr == -Inf] <- 0
log2_expr[is.na(log2_expr)] <- 0

# Output
cat("Writing output to:", output_file, "\n")
write.table(log2_expr, file = output_file, sep = "\t", quote = FALSE, 
            row.names = TRUE, col.names = NA)
