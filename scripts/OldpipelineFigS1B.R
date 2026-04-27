# ---------------------------
# Figure S1B recreation
# ---------------------------

library(readr)
library(dplyr)
library(tibble)
library(FactoMineR)
library(gplots)
library(RColorBrewer)


# input files preparation
expr_matrix <- read_tsv("D:/GroupA_University2026_Project/PCA_Output/2026_UoL_prj_gene_matrix_f.tsv")
metadata <- read_csv("D:/GroupA_University2026_Project/.vs/metadata.csv")


expr_df <- as.data.frame(expr_matrix)
rownames(expr_df) <- expr_df[, 1]          
expr_df <- expr_df[, -1]                   

expr_mat <- as.matrix(expr_df)
mode(expr_mat) <- "numeric"

cat("Total cells before filtering:", ncol(expr_mat), "\n")

# Filter samples using week data instead of tissue
meta_sel <- metadata %>%
  filter(Run %in% colnames(expr_mat)) %>%
  select(Run, stage) %>%
  distinct()

# Extract numeric part from SRR run IDs
meta_sel <- meta_sel %>%
  mutate(
    run_num = as.numeric(gsub("^SRR", "", Run)),
    group = case_when(
      stage == "13 weeks post-conception" ~ "Wpc13",
      run_num >= 2967689 & run_num <= 2967771 ~ "Wpc12 Chip2",
      stage == "12 weeks post-conception" ~ "Wpc12 Chip1",
      TRUE ~ "Other"
    )
  ) %>%
  column_to_rownames("Run")

groups_keep <- c("Wpc12 Chip1", "Wpc12 Chip2", "Wpc13")
meta_sel <- meta_sel[meta_sel$group %in% groups_keep, , drop = FALSE]

samples_keep <- rownames(meta_sel)
expr_mat <- expr_mat[, samples_keep, drop = FALSE]

cat("Samples retained after filtering:", ncol(expr_mat), "\n")
cat("Groups present:", paste(unique(meta_sel$group), collapse = ", "), "\n")

# Variable genes filter

gene_vars <- apply(expr_mat, 1, var)
gene_expr_cells <- rowSums(expr_mat > 1)

variable_genes <- rownames(expr_mat)[
  gene_vars > 0.5 & gene_expr_cells >= 2
]

cat("Variable genes selected:", length(variable_genes), "\n")

if (length(variable_genes) < 10) {
  warning("Very few variable genes – consider adjusting thresholds.")
}

expr_var <- expr_mat[variable_genes, ]


# PCA on filtered samples only
pca_res <- PCA(t(expr_var), scale.unit = TRUE, graph = FALSE)
loadings <- pca_res$var$coord


# Select genes for three sections
n_max <- 50   

# PC1 positive
pc1_pos <- loadings[, 1]
pc1_pos_genes <- names(sort(pc1_pos[pc1_pos > 0.2], decreasing = TRUE))
if (length(pc1_pos_genes) > n_max) pc1_pos_genes <- pc1_pos_genes[1:n_max]

# PC1 negative
pc1_neg <- loadings[, 1]
pc1_neg_genes <- names(sort(pc1_neg[pc1_neg < -0.2], decreasing = FALSE))
if (length(pc1_neg_genes) > n_max) pc1_neg_genes <- pc1_neg_genes[1:n_max]

# PC2 positive
pc2_pos <- loadings[, 2]
pc2_pos_genes <- names(sort(pc2_pos[pc2_pos > 0.2], decreasing = TRUE))
if (length(pc2_pos_genes) > n_max) pc2_pos_genes <- pc2_pos_genes[1:n_max]

cat("\nGenes selected per section:\n")
cat("  PC1 positive:", length(pc1_pos_genes), "\n")
cat("  PC1 negative:", length(pc1_neg_genes), "\n")
cat("  PC2 positive:", length(pc2_pos_genes), "\n")

# Within‑PC hierarchical clustering
cluster_section <- function(gene_list, expr_sub) {
  if (length(gene_list) <= 1) return(gene_list)
  # Compute correlation distance among genes
  gene_cor <- cor(t(expr_sub[gene_list, , drop = FALSE]), method = "pearson")
  gene_dist <- as.dist(1 - gene_cor)
  hc_genes <- hclust(gene_dist, method = "average")
  return(gene_list[hc_genes$order])
}

pc1_pos_ordered <- cluster_section(pc1_pos_genes, expr_mat)
pc1_neg_ordered <- cluster_section(pc1_neg_genes, expr_mat)
pc2_pos_ordered <- cluster_section(pc2_pos_genes, expr_mat)

heat_genes <- c(pc1_pos_ordered, pc1_neg_ordered, pc2_pos_ordered)

cat("Total heatmap genes:", length(heat_genes), "\n")


# Custom colours and labels
meta_ordered <- meta_sel[colnames(expr_mat), , drop = FALSE]
group_colors <- c("Wpc12 Chip1" = "#90EE90",
                  "Wpc12 Chip2" = "#006400",
                  "Wpc13"       = "#6B8E23")
col_side <- group_colors[meta_ordered$group]
col_side <- as.character(col_side)

row_section <- c(
  rep("PC1 pos", length(pc1_pos_ordered)),
  rep("PC1 neg", length(pc1_neg_ordered)),
  rep("PC2 pos", length(pc2_pos_ordered))
)
row_colors <- c("PC1 pos" = "darkgreen", "PC1 neg" = "purple", "PC2 pos" = "darkorange")
row_side <- row_colors[row_section]


# Hierarchical clustering
expr_heat <- expr_mat[heat_genes, ]
cor_mat <- cor(expr_heat, method = "pearson")
cor_dist <- as.dist(1 - cor_mat)
hc_col <- hclust(cor_dist, method = "average")


# Heatmap

hm_colors <- colorRampPalette(c("darkblue", "white", "red"))(100)

row_sep <- cumsum(c(length(pc1_pos_ordered), length(pc1_neg_ordered)))

heatmap.2(expr_heat,
          Colv = as.dendrogram(hc_col),   
          Rowv = FALSE,                    
          dendrogram = "column",           
          col = hm_colors,
          trace = "none",
          density.info = "none",
          key = TRUE,
          keysize = 1.0,
          key.title = "log2(FPKM)",
          margins = c(2, 8),              
          labRow = "",                     
          labCol = "",                     
          ColSideColors = col_side,       
          RowSideColors = row_side,        
          rowsep = row_sep,               
          sepcolor = "black",
          sepwidth = c(0.5, 0.5),
          main = "Figure S1B: Three gene sections with internal clustering")

# Add functional GO annotation text

n_pos <- length(pc1_pos_ordered)   
n_neg <- length(pc1_neg_ordered)   
n_pc2 <- length(pc2_pos_ordered)   

mid_pos <- n_pos / 2
mid_neg <- n_pos + n_neg / 2
mid_pc2 <- n_pos + n_neg + n_pc2 / 2

text_labels <- c(
  "GO: Cell dev., Forebrain dev.",   # PC1 positive
  "GO: Neuron differentiation",     # PC1 negative
  "GO: Cell adhesion, Vesicle transport" # PC2 positive
)

x_text <- 0.15   
y_mid_pos <- 0.68  
y_mid_neg <- 0.35  
y_mid_pc2 <- -0.024 

text(x = x_text, y = y_mid_pos, labels = text_labels[1], 
     srt = 90, xpd = TRUE, col = "white", cex = 0.6, font = 2)
text(x = x_text, y = y_mid_neg, labels = text_labels[2], 
     srt = 90, xpd = TRUE, col = "white", cex = 0.6, font = 2)
text(x = x_text, y = y_mid_pc2, labels = text_labels[3], 
     srt = 90, xpd = TRUE, col = "white", cex = 0.55, font = 2)

# Legend for sample groups
legend("left",
       legend = names(group_colors),
       fill = group_colors,
       border = FALSE,
       bty = "n",
       cex = 0.8,
       title = "Sample Groups")

# Save heatmap
heatmap_plot <- recordPlot()
saveRDS(heatmap_plot, "Figure_S1B_heatmap_plot.rds")

# Write heatmap genes to a text file
output_file <- "heatmap_genes.txt"

file_conn <- file(output_file, "w")

# Write header
writeLines("Genes used in Figure S1B heatmap", file_conn)
writeLines("===================================", file_conn)
writeLines("", file_conn)

writeLines("1. PC1 POSITIVE CORRELATIONS", file_conn)
writeLines(paste0("   Number of genes: ", length(pc1_pos_ordered)), file_conn)
writeLines("   Genes (in clustered order):", file_conn)
for (i in seq_along(pc1_pos_ordered)) {
  writeLines(paste0("   ", i, ". ", pc1_pos_ordered[i]), file_conn)
}
writeLines("", file_conn)

writeLines("2. PC1 NEGATIVE CORRELATIONS (ANTICORRELATED)", file_conn)
writeLines(paste0("   Number of genes: ", length(pc1_neg_ordered)), file_conn)
writeLines("   Genes (in clustered order):", file_conn)
for (i in seq_along(pc1_neg_ordered)) {
  writeLines(paste0("   ", i, ". ", pc1_neg_ordered[i]), file_conn)
}
writeLines("", file_conn)

writeLines("3. PC2 POSITIVE CORRELATIONS", file_conn)
writeLines(paste0("   Number of genes: ", length(pc2_pos_ordered)), file_conn)
writeLines("   Genes (in clustered order):", file_conn)
for (i in seq_along(pc2_pos_ordered)) {
  writeLines(paste0("   ", i, ". ", pc2_pos_ordered[i]), file_conn)
}
writeLines("", file_conn)

writeLines("SUMMARY", file_conn)
writeLines(paste0("Total genes: ", length(heat_genes)), file_conn)
writeLines(paste0(" - PC1 positive: ", length(pc1_pos_ordered)), file_conn)
writeLines(paste0(" - PC1 negative: ", length(pc1_neg_ordered)), file_conn)
writeLines(paste0(" - PC2 positive: ", length(pc2_pos_ordered)), file_conn)

close(file_conn)