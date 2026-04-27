# -----------------------------------------------------------
# Figure S1A PCA reproduction
# -----------------------------------------------------------

library(FactoMineR)
library(ggplot2)

# Load expression matrix
expr_df <- read.delim(
  "2026_UoL_prj_gene_matrix_f.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rownames(expr_df) <- expr_df$gene
expr_df$gene <- NULL
expr_mat <- as.matrix(expr_df)

# Load metadata
meta <- read.csv(
  "metadata.csv",
  stringsAsFactors = FALSE
)

meta_subset <- meta[
  meta$stage %in% c("12 weeks post-conception",
                    "13 weeks post-conception"),
]

# Assign chip information
meta_subset$Chip <- "1"

chip2_range <- paste0("SRR", 2967689:2967771)

meta_subset$Chip[meta_subset$Run %in% chip2_range] <- "2"


# Subset expression matrix
cells_keep <- meta_subset$Run
expr_sub <- expr_mat[, colnames(expr_mat) %in% cells_keep]

meta_subset <- meta_subset[match(colnames(expr_sub), meta_subset$Run), ]


# expressed in > 2 cells
expressed_cells <- rowSums(expr_sub > 0)
keep_expr <- expressed_cells > 2

# variance > 0.5
gene_var <- apply(expr_sub, 1, var, na.rm = TRUE)
keep_var <- gene_var > 0.5
keep_var[is.na(keep_var)] <- FALSE

variable_genes <- keep_expr & keep_var

expr_var <- expr_sub[variable_genes, ]

cat("Genes retained:", sum(variable_genes), "\n")

expr_for_pca <- t(expr_var)

cat("Cells in PCA:", nrow(expr_for_pca), "\n")
cat("Genes in PCA:", ncol(expr_for_pca), "\n")


# PCA using FactoMineR 
pca_res <- PCA(
  expr_for_pca,
  scale.unit = TRUE,
  graph = FALSE
)

# Extract PCA coordinates
pc_df <- as.data.frame(pca_res$ind$coord[, 1:2])
colnames(pc_df) <- c("PC1", "PC2")
pc_df$cell_id <- rownames(pc_df)

pc_df$stage <- meta_subset$stage
pc_df$Chip <- meta_subset$Chip


# Group labels
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


# Custom colours and shapes (original styling)
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


# Variance explained
eig <- pca_res$eig
pc1_var <- round(eig[1, 2], 1)
pc2_var <- round(eig[2, 2], 1)


# Plot PCA 
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
  
  # Full‑width black line
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
# Assign and save plot
pca_plot <- recordPlot()

saveRDS(pca_plot, file = "Oldpipeline_FigS1A.rds")

# Extract sample IDs for Neurons (PC1 < 0) and NPCs (PC1 > 0) and write them out
neurons_ids <- pc_df$cell_id[pc_df$PC1 < 0]
npcs_ids    <- pc_df$cell_id[pc_df$PC1 > 0]


max_len <- max(length(neurons_ids), length(npcs_ids))
output_df <- data.frame(
  Neurons = c(neurons_ids, rep(NA, max_len - length(neurons_ids))),
  NPCs    = c(npcs_ids,    rep(NA, max_len - length(npcs_ids)))
)

write.table(output_df,
            file = "Sample_lists_Neurons_vs_NPCs.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE,
            na = "")
