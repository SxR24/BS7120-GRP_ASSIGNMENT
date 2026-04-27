# -----------------------------------------------------------
# Figure S1A PCA reproduction – Seurat version
# (adapted to use Seurat’s HVGs + PCA instead of FactoMineR)
# -----------------------------------------------------------

library(Seurat)
library(ggplot2)

# Load Seurat object and extract metadata
seurat_obj <- readRDS("seurat_protienpreprocessed.rds")

# Metadata with Run column
meta <- seurat_obj@meta.data
meta$Run <- rownames(meta)

# Stage & chip filtering
meta_subset <- meta[
  meta$stage %in% c("12 weeks post-conception",
                    "13 weeks post-conception"),
]

meta_subset$Chip <- "1"
chip2_range <- paste0("SRR", 2967689:2967771)
meta_subset$Chip[meta_subset$Run %in% chip2_range] <- "2"

cells_keep <- meta_subset$Run

# Subset Seurat object to these cells only
obj_sub <- subset(seurat_obj, cells = cells_keep)

# Highly variable genes (VST, n = 3000) and PCA
obj_sub <- FindVariableFeatures(obj_sub,
                                selection.method = "vst",
                                nfeatures = 3000)
obj_sub <- ScaleData(obj_sub, features = VariableFeatures(obj_sub))
obj_sub <- RunPCA(obj_sub, features = VariableFeatures(obj_sub))

# Extract PCA coordinates and variance explained
pca_coords <- Embeddings(obj_sub, "pca")[, 1:2]
variance_explained <- (obj_sub[["pca"]]@stdev^2) /
  sum(obj_sub[["pca"]]@stdev^2) * 100
pc1_var <- round(variance_explained[1], 1)
pc2_var <- round(variance_explained[2], 1)

cat("Genes retained (HVGs):", length(VariableFeatures(obj_sub)), "\n")
cat("Cells in PCA:", nrow(pca_coords), "\n")
cat("Genes in PCA:", length(VariableFeatures(obj_sub)), "\n")

# Build plotting dataframe
pc_df <- data.frame(
  PC1     = pca_coords[, 1],
  PC2     = pca_coords[, 2],
  cell_id = rownames(pca_coords),
  stringsAsFactors = FALSE
)

pc_df$stage <- meta_subset$stage
pc_df$Chip  <- meta_subset$Chip

pc_df$Group <- NA
pc_df$Group[pc_df$stage == "13 weeks post-conception"] <- "13 wpc"
pc_df$Group[pc_df$stage == "12 weeks post-conception" &
              pc_df$Chip == "2"] <- "12 wpc chip 2"
pc_df$Group[pc_df$stage == "12 weeks post-conception" &
              pc_df$Chip == "1"] <- "12 wpc chip 1"
pc_df$Group <- factor(
  pc_df$Group,
  levels = c("13 wpc",
             "12 wpc chip 2",
             "12 wpc chip 1")
)

# Plotting

my_colours <- c(
  "13 wpc" = "#6B8E23",
  "12 wpc chip 2" = "#006400",
  "12 wpc chip 1" = "#90EE90"
)
my_shapes <- c(
  "13 wpc" = 16,
  "12 wpc chip 2" = 16,
  "12 wpc chip 1" = 15
)

y_min <- floor(min(pc_df$PC2))
y_max <- ceiling(max(pc_df$PC2))
y_range <- y_max - y_min
y_line <- y_max + y_range * 0.07

ggplot(pc_df, aes(x = PC1, y = PC2, colour = Group, shape = Group)) +
  geom_point(size = 3) +
  scale_colour_manual(values = my_colours) +
  scale_shape_manual(values = my_shapes) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    colour = "",
    shape = ""
  ) +
  scale_x_continuous(breaks = seq(-20, 100, by = 20)) +
  scale_y_continuous(breaks = seq(-20, 100, by = 20)) +
  theme_light(base_size = 14) +
  theme(
    legend.position = "right",
    plot.margin = margin(t = 35, r = 10, b = 10, l = 10)
  ) +
  # Full‑width black line with arrows
  geom_segment(
    aes(x = -Inf, xend = Inf, y = y_line, yend = y_line),
    colour = "black", size = 1.2,
    arrow = arrow(ends = "both", length = unit(0.18, "cm"), type = "closed"),
    inherit.aes = FALSE
  ) +
  # Black dot at zero
  geom_point(
    aes(x = 0, y = y_line),
    colour = "black", size = 3, shape = 16,
    inherit.aes = FALSE
  ) +
  # Neuron / NPC labels
  geom_text(
    aes(x = -5, y = y_line, label = "Neurons"),
    colour = "black", size = 5,
    hjust = 1, vjust = -1.0,
    inherit.aes = FALSE
  ) +
  geom_text(
    aes(x = 5, y = y_line, label = "NPCs"),
    colour = "black", size = 5,
    hjust = 0, vjust = -1.0,
    inherit.aes = FALSE
  ) +
  coord_cartesian(ylim = c(y_min, y_max), clip = "off")

# Save plot
pca_plot <- recordPlot()
saveRDS(pca_plot, file = "Oldpipeline_FigS1A.rds")

# Export cell IDs
neurons_ids <- pc_df$cell_id[pc_df$PC1 < 0]
npcs_ids    <- pc_df$cell_id[pc_df$PC1 > 0]

max_len <- max(length(neurons_ids), length(npcs_ids))
output_df <- data.frame(
  Neurons = c(neurons_ids, rep(NA, max_len - length(neurons_ids))),
  NPCs    = c(npcs_ids,    rep(NA, max_len - length(npcs_ids)))
)

write.table(output_df,
            file = "protien_Sample_lists_Neurons_vs_NPCs.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE,
            na = "")