# =============================================================================
# Monocle 3 with ICA
# =============================================================================

library(monocle3)
library(dplyr)
library(fastICA)
library(ggplot2)

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

# Read zone assignments
assign_file <- "cell_assignments.txt"
assign_df <- read.table(assign_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
stopifnot(all(c("Sample", "Assigned_Zone") %in% colnames(assign_df)))
assign_df <- assign_df %>% filter(Sample %in% colnames(expr_mat))

cell_meta <- data.frame(cell = colnames(expr_mat),
                        row.names = colnames(expr_mat))
cell_meta$zone <- assign_df$Assigned_Zone[match(rownames(cell_meta), assign_df$Sample)]

cds <- new_cell_data_set(expr_mat,
                         cell_metadata = cell_meta,
                         gene_metadata = data.frame(gene_short_name = rownames(expr_mat),
                                                    row.names = rownames(expr_mat)))

# Dummy size factors
cds@colData$Size_Factor <- rep(1, ncol(expr_mat))

# Run ICA and store as UMAP
set.seed(12345)
ica_res <- fastICA(t(expr_mat), n.comp = 2)
ica_coords <- ica_res$S
rownames(ica_coords) <- colnames(expr_mat)
colnames(ica_coords) <- c("IC1", "IC2")

reducedDim(cds, "UMAP") <- ica_coords

# Cluster cells on ICA coordinates
cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = 0.007)

#Learn graph (MST) in ICA space
cds <- learn_graph(cds,
                   use_partition = FALSE,
                   close_loop = FALSE,
                   learn_graph_control = list(minimal_branch_len = 1))

zone_colors <- c(
  "VZ"  = "#FFBF20",
  "iSVZ"= "#F77800",
  "oSVZ"= "#CD2624",
  "CP"  = "#7D007E"
)

# Plot
plot_cells(cds,
           color_cells_by = "zone",
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 1.5,
           cell_size = 1.5) +
  scale_color_manual(values = zone_colors) +
  theme(legend.position = "right") +
  ggtitle("Fetal Neocortex Lineage Tree (ICA + MST) – Zones")

# =============================================================================
# Replicate Fig2B
# =============================================================================

# Define the genes of interest
genes_of_interest <- c("SOX2", "EOMES", "MYT1L")

available_genes <- genes_of_interest[genes_of_interest %in% rownames(expr_mat)]
if (length(available_genes) < length(genes_of_interest)) {
  missing_genes <- setdiff(genes_of_interest, available_genes)
  warning("The following genes are not in the expression matrix and will be skipped: ",
          paste(missing_genes, collapse = ", "))
}

# Loop over available genes and produce a plot for each
for (gene in available_genes) {
  
  colData(cds)[[paste0("raw_", gene)]] <- as.numeric(expr_mat[gene, colnames(cds)])
  
  p <- plot_cells(cds,
                  color_cells_by = paste0("raw_", gene),
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  graph_label_size = 1.5,
                  cell_size = 1.5)
  
  p <- p +
    scale_color_gradient(low = "darkblue", high = "orange", name = "Raw expression") +
    theme(legend.position = "right") +
    ggtitle(paste("Fetal Neocortex Lineage Tree (ICA + MST) –", gene, "Raw Expression"))
  
  print(p)
  
}