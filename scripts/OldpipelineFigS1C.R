# -------------------------------------
# Corrected Replication of Figure S1C
# -------------------------------------

library(readr)
library(dplyr)
library(tibble)
library(gplots)
library(RColorBrewer)
library(stats)

expr_matrix <- read_tsv("D:/GroupA_University2026_Project/PCA_Output/2026_UoL_prj_gene_matrix_f.tsv")
metadata <- read_csv("D:/GroupA_University2026_Project/.vs/metadata.csv")

expr_df <- as.data.frame(expr_matrix)
rownames(expr_df) <- expr_df[, 1]
expr_df <- expr_df[, -1]
expr_mat <- as.matrix(expr_df)
mode(expr_mat) <- "numeric"

meta_clean <- metadata %>%
  filter(stage %in% c("12 weeks post-conception", "13 weeks post-conception")) %>%
  select(Run, stage) %>%
  distinct() %>%
  mutate(
    chip = case_when(
      stage == "13 weeks post-conception" ~ "Wpc13",
      Run >= "SRR2967689" & Run <= "SRR2967771" ~ "Wpc12 Chip2",
      stage == "12 weeks post-conception" ~ "Wpc12 Chip1",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(chip))

keep_samples <- intersect(meta_clean$Run, colnames(expr_mat))
expr_filt <- expr_mat[, keep_samples]

hk_genes <- c("ACTB", "GAPDH")
hk_present <- hk_genes[hk_genes %in% rownames(expr_filt)]
if (length(hk_present) == 0) stop("ACTB or GAPDH not found.")
keep_cells <- colnames(expr_filt)[colSums(expr_filt[hk_present, , drop = FALSE] > 0) > 0]
expr_filt <- expr_filt[, keep_cells]
meta_clean <- meta_clean %>% filter(Run %in% keep_cells)

marker_genes <- c("PAX6", "SOX2", "VIM", "MYT1L", "TBR1",
                  "ERBB4", "DLX1", "DLX2", "GAD1", "GLI3",
                  "NEUROD6", "DLX5", "DLX6")
present_genes <- marker_genes[marker_genes %in% rownames(expr_filt)]
expr_markers <- expr_filt[present_genes, , drop = FALSE]

gene_vars <- apply(expr_markers, 1, var)
zero_var_genes <- names(gene_vars[gene_vars == 0 | is.na(gene_vars)])
if (length(zero_var_genes) > 0) {
  expr_markers <- expr_markers[!rownames(expr_markers) %in% zero_var_genes, , drop = FALSE]
}

ann_colors <- list(
  chip = c("Wpc12 Chip1" = "#90EE90",
           "Wpc12 Chip2" = "#006400",
           "Wpc13" = "#6B8E23")
)
meta_clean <- meta_clean %>% column_to_rownames("Run")
meta_clean <- meta_clean[colnames(expr_markers), , drop = FALSE]
col_side_colors <- cbind(Chip = ann_colors$chip[meta_clean$chip])

# Cell clustering
cor_mat <- cor(expr_markers, method = "pearson")
cor_mat[is.na(cor_mat)] <- 0
cor_dist <- as.dist(1 - cor_mat)
hc_cells <- hclust(cor_dist, method = "average")

# Gene clustering
gene_cor_mat <- cor(t(expr_markers), method = "pearson")
gene_cor_mat[is.na(gene_cor_mat)] <- 0
gene_cor_dist <- as.dist(1 - gene_cor_mat)
hc_genes <- hclust(gene_cor_dist, method = "average")

# Heatmap
hm_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

heatmap.2(expr_markers,
          Colv = as.dendrogram(hc_cells),
          Rowv = as.dendrogram(hc_genes),    # preserves exact leaf order
          dendrogram = "column",             # hides row dendrogram
          col = hm_colors,
          trace = "none",
          density.info = "none",
          key = TRUE,
          keysize = 1.0,
          key.title = "log2(FPKM)",
          margins = c(12, 10),
          cexRow = 0.9,
          labCol = "",                       # hide column labels
          ColSideColors = col_side_colors,
          main = "Figure S1C: Interneuron marker expression",
          xlab = "Cells",
          ylab = "Marker Genes")

legend("topright",
       legend = names(ann_colors$chip),
       fill = ann_colors$chip,
       border = FALSE, bty = "n", cex = 0.8)

# Interneuron identification
int_markers_sub <- c("DLX1", "DLX2", "DLX5", "DLX6", "GAD1", "ERBB4")
int_present <- int_markers_sub[int_markers_sub %in% rownames(expr_markers)]
if (length(int_present) > 0) {
  int_score <- colMeans(expr_markers[int_present, , drop = FALSE])
  top_int_cells <- names(sort(int_score, decreasing = TRUE))[1:5]
  cat("\nPredicted interneurons (top 5 scoring cells):\n")
  print(top_int_cells)
}

heatmap_plot <- recordPlot()
saveRDS(heatmap_plot, "D:/GroupA_University2026_Project/OUTPUTS/Figure_S1C_heatmap_plot.rds")