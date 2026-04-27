# =============================================================================
# Corrected Replication of Figure S1C (adapated for seurat)
# =============================================================================

library(Seurat)
library(readr)
library(dplyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(stats)

# Input
seurat_preprocessing <- readRDS("seurat_preprocessing")
expr_mat_full <- as.matrix(GetAssayData(seurat_preprocessed.rds, layer = "data"))

# Extract metadata
meta_data <- seurat_protienpreprocessing@meta.data
meta_clean <- data.frame(
  Run   = colnames(seurat_protienpreprocessing),   
  stage = meta_data$stage,                        
  stringsAsFactors = FALSE
) %>%
  filter(stage %in% c("12 weeks post-conception", "13 weeks post-conception")) %>%
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

common_cells <- intersect(meta_clean$Run, colnames(expr_mat_full))
expr_filt <- expr_mat_full[, common_cells, drop = FALSE]
meta_clean <- meta_clean %>% filter(Run %in% common_cells)

# Marker genes assigned

marker_genes <- c("PAX6", "SOX2", "VIM", "MYT1L", "TBR1",
                  "ERBB4", "DLX1", "DLX2", "GAD1", "GLI3",
                  "NEUROD6", "DLX5", "DLX6")
present_genes <- marker_genes[marker_genes %in% rownames(expr_filt)]
expr_markers <- expr_filt[present_genes, , drop = FALSE]

# Prepare visual affects

chip_colors <- c("Wpc12 Chip1" = "#90EE90",
                 "Wpc12 Chip2" = "#006400",
                 "Wpc13"      = "#6B8E23")

meta_clean <- meta_clean %>% column_to_rownames("Run")
meta_clean <- meta_clean[colnames(expr_markers), , drop = FALSE]
col_anno <- data.frame(Chip = meta_clean$chip, row.names = rownames(meta_clean))

ann_colors <- list(Chip = chip_colors)

# Pheatmap

my_colours <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

pheatmap_obj <- pheatmap(expr_markers,
                         scale = "none",               # already log‑normalised
                         clustering_method = "ward.D2",
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         show_rownames = TRUE,
                         show_colnames = FALSE,        # too many cells to label
                         color = my_colours,
                         main = "Figure S1C: Interneuron marker expression",
                         fontsize_row = 8,
                         border_color = NA,
                         treeheight_row = 50,
                         treeheight_col = 50)

saveRDS(pheatmap_obj, "D:/GroupA_University2026_Project/OUTPUTS/AltFigure_S1C_heatmap_protien_pheatmap.rds")

# Interneuron identification

int_markers_sub <- c("DLX1", "DLX2", "DLX5", "DLX6", "GAD1", "ERBB4")
int_present <- int_markers_sub[int_markers_sub %in% rownames(expr_markers)]
if (length(int_present) > 0) {
  int_score <- colMeans(expr_markers[int_present, , drop = FALSE])
  top_int_cells <- names(sort(int_score, decreasing = TRUE))[1:5]
  cat("\nPredicted interneurons (top 5 scoring cells):\n")
  print(top_int_cells)
}