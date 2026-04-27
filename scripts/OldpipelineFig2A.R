# =============================================================================
# Monocle 3 with ICA
# =============================================================================

library(monocle3)
library(dplyr)
library(fastICA)

# Load and filter data
mat_file <- "2026_UoL_prj_gene_matrix_f.tsv"
expr_mat <- as.matrix(read.table(mat_file, header = TRUE, row.names = 1, 
                                 sep = "\t", check.names = FALSE))

gene_list_file <- "Monoclegenelist.csv"
gene_df <- read.csv(gene_list_file, header = TRUE, stringsAsFactors = FALSE)
genes_to_keep <- intersect(gene_df[[1]], rownames(expr_mat))
expr_mat <- expr_mat[genes_to_keep, ]

metadata_file <- "metadata.csv"
meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
fetal_cells <- meta %>% filter(tissue == "Fetal neocortex") %>% pull(Run)
exclude_runs <- c("SRR2967628", "SRR2967645", "SRR2967696", "SRR2967704", "SRR2967725")
fetal_cells <- setdiff(fetal_cells, exclude_runs)
cells_to_keep <- intersect(fetal_cells, colnames(expr_mat))
expr_mat <- expr_mat[, cells_to_keep]

assign_file <- "oldsubpopulations.csv"
assign_df <- read.csv(assign_file, stringsAsFactors = FALSE)
colnames(assign_df)[colnames(assign_df) == "sample_id"] <- "sample_id"
colnames(assign_df)[colnames(assign_df) == "group"] <- "group"
assign_df <- assign_df %>% filter(sample_id %in% colnames(expr_mat))

cell_meta <- data.frame(cell = colnames(expr_mat), row.names = colnames(expr_mat))
cell_meta$subpopulation <- assign_df$group[match(rownames(cell_meta), assign_df$sample_id)]

cds <- new_cell_data_set(expr_mat,
                         cell_metadata = cell_meta,
                         gene_metadata = data.frame(gene_short_name = rownames(expr_mat),
                                                    row.names = rownames(expr_mat)))

# Dummy size factors
cds@colData$Size_Factor <- rep(1, ncol(expr_mat))

# Run ICA
set.seed(12345)
ica_res <- fastICA(t(expr_mat), n.comp = 2)
ica_coords <- ica_res$S
rownames(ica_coords) <- colnames(expr_mat)
colnames(ica_coords) <- c("IC1", "IC2")

# Store ICA as the "UMAP" reduction
reducedDim(cds, "UMAP") <- ica_coords

# Cluster cells on ICA coordinates
cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = 0.007)

# Learn graph (MST) in ICA space
cds <- learn_graph(cds,
                   use_partition = FALSE,      
                   close_loop = FALSE,         
                   learn_graph_control = list(minimal_branch_len = 1))

group_cols <- c(
  "AP2" = "#95F297",
  "AP1" = "#2C8334",
  "BP1" = "#15B5B5",
  "BP2" = "#9DFFFF",
  "N1"  = "#8EC9F5",
  "N2"  = "#4589EF",
  "N3"  = "#0545BB"
)

# Plot
plot_cells(cds,
           color_cells_by = "subpopulation",
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 1.5,
           cell_size = 1.5) +
  scale_color_manual(values = group_cols) +
  theme(legend.position = "right") +
  ggtitle("Fetal Neocortex Lineage Tree (ICA + MST)")