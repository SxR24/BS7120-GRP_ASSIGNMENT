# =============================================================================
# Figure 1B 
# =============================================================================

library(readr)
library(dplyr)
library(tibble)
library(gplots)
library(stats)

# Inputs
sc_expr_file    <- "D:/GroupA_University2026_Project/PCA_Output/2026_UoL_prj_gene_matrix_f.tsv"
metadata_file   <- "D:/GroupA_University2026_Project/.vs/metadata.csv"
bulk_zones_file <- "D:/GroupA_University2026_Project/BulkRNA/Old_GSE38805_human_FPKM_matrix_log2_symbols.txt"
bulk_celltypes_file <- "D:/GroupA_University2026_Project/BulkRNA/Old_GSE65000_hsa_fpkm_matrix_log2.txt"

interneurons_to_remove <- c("SRR2967704", "SRR2967725", "SRR2967696", 
                            "SRR2967645", "SRR2967628")

# Load single-cell data and metadata
sc_raw <- read_tsv(sc_expr_file)
sc_df <- as.data.frame(sc_raw)
rownames(sc_df) <- sc_df[, 1]
sc_df <- sc_df[, -1]
sc_mat <- as.matrix(sc_df)
mode(sc_mat) <- "numeric"

metadata <- read_csv(metadata_file)

# keep only week 12 and week 13 samples
meta_filt <- metadata %>%
  filter(stage %in% c("12 weeks post-conception", "13 weeks post-conception")) %>%
  select(Run, stage) %>%
  distinct()

# Keep only samples that exist in the expression matrix
keep_samples <- intersect(meta_filt$Run, colnames(sc_mat))
sc_mat <- sc_mat[, keep_samples, drop = FALSE]

cat(sprintf("Samples retained after early filter (week 12 & 13 only): %d\n", ncol(sc_mat)))

# Remove the five interneurons
sc_mat <- sc_mat[, !colnames(sc_mat) %in% interneurons_to_remove, drop = FALSE]
cat(sprintf("Cells remaining after interneuron removal: %d\n", ncol(sc_mat)))

# Load and process bulk zone data (GSE38805)
bulk_zones <- read_tsv(bulk_zones_file)
bulk_zones_df <- as.data.frame(bulk_zones)
rownames(bulk_zones_df) <- bulk_zones_df[, 1]
bulk_zones_df <- bulk_zones_df[, -1]
bulk_zones_mat <- as.matrix(bulk_zones_df)
mode(bulk_zones_mat) <- "numeric"

# Extract only week 13 columns
w13_cols <- grep("^w13_", colnames(bulk_zones_mat), value = TRUE)
cat("Week 13 zone columns found:\n")
print(w13_cols)

get_zone_mean <- function(zone) {
  cols <- grep(paste0("w13_", zone, "$"), w13_cols, value = TRUE)
  if (length(cols) == 0) stop("No column for zone ", zone)
  rowMeans(bulk_zones_mat[, cols, drop = FALSE])
}

zones <- c("VZ", "ISVZ", "OSVZ", "CP")
bulk_zones_avg <- do.call(cbind, lapply(zones, get_zone_mean))
colnames(bulk_zones_avg) <- c("VZ", "iSVZ", "oSVZ", "CP")
rownames(bulk_zones_avg) <- rownames(bulk_zones_mat)
bulk_zones_avg <- bulk_zones_avg[complete.cases(bulk_zones_avg), , drop = FALSE]

# Load and process bulk cell type data (GSE65000)
bulk_celltypes <- read_tsv(bulk_celltypes_file)
bulk_celltypes_df <- as.data.frame(bulk_celltypes)
rownames(bulk_celltypes_df) <- bulk_celltypes_df[, 1]
bulk_celltypes_df <- bulk_celltypes_df[, -1]
bulk_celltypes_mat <- as.matrix(bulk_celltypes_df)
mode(bulk_celltypes_mat) <- "numeric"

get_type_mean <- function(type) {
  cols <- grep(paste0("Hsa_", type, "_"), colnames(bulk_celltypes_mat), value = TRUE)
  if (length(cols) == 0) stop("No columns for type ", type)
  rowMeans(bulk_celltypes_mat[, cols, drop = FALSE])
}

celltypes <- c("aRG", "bRG", "N")
bulk_celltypes_avg <- do.call(cbind, lapply(celltypes, get_type_mean))
colnames(bulk_celltypes_avg) <- c("aRG", "bRG", "neuron")
rownames(bulk_celltypes_avg) <- rownames(bulk_celltypes_mat)
bulk_celltypes_avg <- bulk_celltypes_avg[complete.cases(bulk_celltypes_avg), , drop = FALSE]

# match genes between single-cell and bulk
common_genes_zones <- intersect(rownames(sc_mat), rownames(bulk_zones_avg))
common_genes_celltypes <- intersect(rownames(sc_mat), rownames(bulk_celltypes_avg))

cat(sprintf("Common genes with zones: %d\n", length(common_genes_zones)))
cat(sprintf("Common genes with cell types: %d\n", length(common_genes_celltypes)))

sc_zones <- sc_mat[common_genes_zones, , drop = FALSE]
bulk_zones_sub <- bulk_zones_avg[common_genes_zones, , drop = FALSE]

sc_celltypes <- sc_mat[common_genes_celltypes, , drop = FALSE]
bulk_celltypes_sub <- bulk_celltypes_avg[common_genes_celltypes, , drop = FALSE]

# Spearman correlation per single cell
correlate_cell <- function(cell_expr, bulk_mat) {
  apply(bulk_mat, 2, function(bulk_col) {
    cor(cell_expr, bulk_col, method = "spearman", use = "complete.obs")
  })
}

cor_zones <- t(apply(sc_zones, 2, correlate_cell, bulk_mat = bulk_zones_sub))
cor_celltypes <- t(apply(sc_celltypes, 2, correlate_cell, bulk_mat = bulk_celltypes_sub))

zscore_norm <- function(mat) t(scale(t(mat)))

cor_zones_scaled <- zscore_norm(cor_zones)
cor_celltypes_scaled <- zscore_norm(cor_celltypes)

# Hierarchical clustering
cell_dist_zones <- as.dist(1 - cor(t(cor_zones_scaled), method = "pearson"))
hc_zones <- hclust(cell_dist_zones, method = "average")

cell_dist_celltypes <- as.dist(1 - cor(t(cor_celltypes_scaled), method = "pearson"))
hc_celltypes <- hclust(cell_dist_celltypes, method = "average")

# Metadata for annotation
meta_plot <- meta_filt %>%
  filter(Run %in% colnames(sc_mat)) %>%
  mutate(
    chip = case_when(
      stage == "13 weeks post-conception" ~ "Wpc13",
      Run >= "SRR2967689" & Run <= "SRR2967771" ~ "Wpc12 Chip2",
      stage == "12 weeks post-conception" ~ "Wpc12 Chip1",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(chip)) %>%
  column_to_rownames("Run")

meta_plot <- meta_plot[colnames(sc_mat), , drop = FALSE]

ann_colors <- list(
  chip = c("Wpc12 Chip1" = "#90EE90",
           "Wpc12 Chip2" = "#006400",
           "Wpc13"     = "#6B8E23")
)

row_side_colors <- cbind(Chip = ann_colors$chip[meta_plot$chip])

# 11. Heatmaps
hm_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Zones heatmap
heatmap.2(cor_zones_scaled,
          Colv = NA,
          Rowv = as.dendrogram(hc_zones),
          dendrogram = "row",
          col = hm_colors,
          trace = "none",
          density.info = "none",
          key = TRUE,
          keysize = 1.0,
          key.title = "Z-score",
          margins = c(8, 10),
          cexRow = 0.6,
          cexCol = 1.0,
          labRow = "",
          labCol = colnames(cor_zones_scaled),
          RowSideColors = row_side_colors,
          main = "Correlation with microdissected zones",
          xlab = "Cortical Zones")

legend("topright",
       legend = names(ann_colors$chip),
       fill = ann_colors$chip,
       border = FALSE, bty = "n", cex = 1.2)

zones_plot <- recordPlot()
saveRDS(zones_plot, file = "D:/GroupA_University2026_Project/OUTPUTS/Figure_1B_heatmapa.rds")

# Cell types heatmap
heatmap.2(cor_celltypes_scaled,
          Colv = NA,
          Rowv = as.dendrogram(hc_celltypes),
          dendrogram = "row",
          col = hm_colors,
          trace = "none",
          density.info = "none",
          key = TRUE,
          keysize = 1.0,
          key.title = "Z-score",
          margins = c(8, 10),
          cexRow = 0.6,
          cexCol = 1.0,
          labRow = "",
          labCol = colnames(cor_celltypes_scaled),
          RowSideColors = row_side_colors,
          main = "Correlation with FACS-purified cell types",
          xlab = "Cell Types")

legend("topright",
       legend = names(ann_colors$chip),
       fill = ann_colors$chip,
       border = FALSE, bty = "n", cex = 1.2)

zones_plot <- recordPlot()
saveRDS(zones_plot, file = "D:/GroupA_University2026_Project/OUTPUTS/Figure_1B_heatmapb.rds")

# Save results
results_fig1b <- list(
  cor_zones = cor_zones,
  cor_zones_scaled = cor_zones_scaled,
  cor_celltypes = cor_celltypes,
  cor_celltypes_scaled = cor_celltypes_scaled,
  hc_zones = hc_zones,
  hc_celltypes = hc_celltypes,
  metadata = meta_plot,
  sc_mat_filtered = sc_mat
)
saveRDS(results_fig1b, file = "D:/GroupA_University2026_Project/OUTPUTS/Figure_1B_results.rds")

cat("\nFigure 1B completed. Results saved to Figure_1B_results.rds\n")

# Create assignment table: sample for assigned_celltype and assigned_zone

zone_assignment <- colnames(cor_zones)[max.col(cor_zones, ties.method = "first")]
names(zone_assignment) <- rownames(cor_zones)

celltype_assignment <- colnames(cor_celltypes)[max.col(cor_celltypes, ties.method = "first")]
names(celltype_assignment) <- rownames(cor_celltypes)

common_cells <- intersect(names(zone_assignment), names(celltype_assignment))

# Create a data frame
assignment_df <- data.frame(
  Sample = common_cells,
  Assigned_Zone = zone_assignment[common_cells],
  Assigned_CellType = celltype_assignment[common_cells],
  stringsAsFactors = FALSE
)

assignment_df <- assignment_df %>%
  left_join(meta_plot %>% rownames_to_column("Sample") %>% select(Sample, chip), by = "Sample")

write.table(assignment_df,
            file = "D:/GroupA_University2026_Project/OUTPUTS/cell_assignments.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

cat("Assignment table saved to cell_assignments.txt\n")

print(head(assignment_df))