# ---------------------------------------------------------------
# Reconstruct Figure 3D 
# Gene‑centric PC selection + phylogenetic tree + ROC merging
# ---------------------------------------------------------------

library(Seurat)          
library(FactoMineR)      
library(dplyr)
library(ggplot2)
library(Matrix)
library(scales)
library(tidyr)
library(patchwork)


# Load and filter data
mat <- read.table("D:/GroupA_University2026_Project/PCA_Output/2026_UoL_prj_gene_matrix_f.tsv",
                  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

meta <- read.csv("D:/GroupA_University2026_Project/.vs/metadata.csv",
                 stringsAsFactors = FALSE)

region_meta <- read.csv("D:/GroupA_University2026_Project/OUTPUTS/Organoidlabelmetadata.csv",
                        stringsAsFactors = FALSE)

organoid_days <- c("33 days", "35 days", "37 days", "41 days", "65 days", "53 days", "58 days")
organoid_meta <- meta[meta$stage %in% organoid_days, ]

cells_to_keep <- intersect(colnames(mat), organoid_meta$Run)
mat_org <- as.matrix(mat[, cells_to_keep])

rownames(organoid_meta) <- organoid_meta$Run
organoid_meta <- organoid_meta[colnames(mat_org), ]

organoid_meta$region <- NA
region_meta_sub <- region_meta[region_meta$sample %in% rownames(organoid_meta), ]
organoid_meta[region_meta_sub$sample, "region"] <- region_meta_sub$group

# Create Seurat object
rna_assay <- CreateAssayObject(counts = mat_org)
rna_assay <- SetAssayData(rna_assay, layer = "data", new.data = mat_org)

seurat_obj <- CreateSeuratObject(rna_assay, meta.data = organoid_meta)
DefaultAssay(seurat_obj) <- "RNA"

#Identify variable genes
expr_data <- LayerData(seurat_obj, layer = "data")
genes_expressed <- rowSums(expr_data > 0) > 3
gene_vars <- apply(expr_data[genes_expressed, ], 1, var)
var_genes <- names(gene_vars[gene_vars > 1])

VariableFeatures(seurat_obj) <- var_genes
message("Selected ", length(var_genes), " variable genes.")

# PCA using FactoMineR
expr_for_pca <- t(expr_data[var_genes, ])   

pca_facto <- PCA(expr_for_pca, scale.unit = TRUE, ncp = 50, graph = FALSE)

# Store PCA results in Seurat object
seurat_obj[["pca"]] <- CreateDimReducObject(
  embeddings = pca_facto$ind$coord,
  loadings  = pca_facto$var$coord,
  stdev     = pca_facto$eig[, 1],
  key       = "PC_",
  assay     = "RNA"
)

# Gene‑level permutation test
set.seed(123)
n_permutations <- 200
n_genes <- length(var_genes)
n_shuffle <- round(n_genes * 0.01)   
max_pcs_to_test <- 20                 

message("Starting gene‑level permutation test (", n_permutations, " permutations)...")

perm_sds <- matrix(NA, nrow = n_permutations, ncol = max_pcs_to_test)
perm_loadings_list <- vector("list", n_permutations)

for (i in seq_len(n_permutations)) {
  shuffle_idx <- sample(1:ncol(expr_for_pca), n_shuffle)
  perm_expr <- expr_for_pca
  perm_expr[, shuffle_idx] <- perm_expr[, sample(shuffle_idx)]
  
  perm_pca <- prcomp(perm_expr, center = TRUE, scale. = TRUE, rank. = max_pcs_to_test)
  
  perm_loadings_list[[i]] <- perm_pca$rotation      
  perm_sds[i, ] <- perm_pca$sdev[1:max_pcs_to_test] 
  
  if (i %% 20 == 0) message("  Permutation ", i, " complete")
}

obs_loadings <- pca_facto$var$coord[, 1:max_pcs_to_test]
observed_sds <- pca_facto$eig[1:max_pcs_to_test, 1]

pval_matrix <- matrix(NA, nrow = n_genes, ncol = max_pcs_to_test,
                      dimnames = list(var_genes, paste0("PC", 1:max_pcs_to_test)))
for (pc in 1:max_pcs_to_test) {
  for (g in seq_len(n_genes)) {
    obs <- abs(obs_loadings[g, pc])
    perm_vals <- sapply(perm_loadings_list, function(x) abs(x[g, pc]))
    pval_matrix[g, pc] <- mean(perm_vals >= obs)
  }
}

# Gene‑centric PC selection
genes_extreme_per_pc <- apply(pval_matrix, 2, function(x) sum(x < 1e-20))
significant_pcs <- which(genes_extreme_per_pc > 0)

if (length(significant_pcs) == 0) {
  warning("No PC contains a gene with P < 1e‑20. Falling back to first 3 PCs.")
  significant_pcs <- 1:3
} else {
  message("PCs with at least one gene P < 1e‑20: ", paste(significant_pcs, collapse = ", "))
}

# Select genes from significant PCs
selected_genes <- c()
for (pc in significant_pcs) {
  gene_pvals <- pval_matrix[, pc]
  sig_genes <- names(gene_pvals)[gene_pvals < 1e-3]
  loadings_abs <- abs(obs_loadings[sig_genes, pc])
  top_genes <- names(sort(loadings_abs, decreasing = TRUE))[1:min(50, length(sig_genes))]
  selected_genes <- union(selected_genes, top_genes)
}
message("Selected ", length(selected_genes), " genes for downstream analysis.")

# PCA on selected genes
seurat_obj <- ScaleData(seurat_obj, features = selected_genes,
                        do.scale = FALSE, do.center = FALSE)
seurat_obj <- RunPCA(seurat_obj, features = selected_genes, npcs = min(50, length(selected_genes)))

# t‑SNE with perplexity = 5
seurat_obj <- RunTSNE(seurat_obj, dims = 1:20, perplexity = 5,
                      check_duplicates = FALSE, seed.use = 42)

# Over‑clustering using PCA
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 2.0)  
message("Initial clusters: ", length(levels(seurat_obj)))

# Build phylogenetic tree 
seurat_obj <- BuildClusterTree(seurat_obj, dims = 1:20)

# visualise the tree
PlotClusterTree(seurat_obj)

# ROC‑based branch collapsing
get_max_roc_auc <- function(obj, id1, id2) {
  markers <- FindMarkers(obj, ident.1 = id1, ident.2 = id2,
                         test.use = "roc", logfc.threshold = 0,
                         min.pct = 0, verbose = FALSE)
  if (nrow(markers) == 0) return(0)
  max(markers$myAUC)
}

tree <- Tool(seurat_obj, "BuildClusterTree")
if (is.null(tree)) stop("Phylogenetic tree not found. Did BuildClusterTree run successfully?")

edges <- tree$edge
parent_nodes <- unique(edges[, 1])
sibling_pairs <- list()

for (p in parent_nodes) {
  children <- edges[edges[, 1] == p, 2]
  if (length(children) == 2 && all(children <= length(tree$tip.label))) {
    cl1 <- tree$tip.label[children[1]]
    cl2 <- tree$tip.label[children[2]]
    sibling_pairs[[length(sibling_pairs) + 1]] <- c(cl1, cl2)
  }
}

cat("\n========== SIBLING CLUSTER PAIRS ==========\n")
for (pair in sibling_pairs) {
  auc <- get_max_roc_auc(seurat_obj, pair[1], pair[2])
  cat(sprintf("Clusters %s and %s  -->  Max ROC AUC = %.3f\n", pair[1], pair[2], auc))
}
cat("===========================================\n\n")

# Merge sibling pairs with Low AUC
auto_merge_threshold <- 0.7 

for (pair in sibling_pairs) {
  if (!all(pair %in% levels(seurat_obj))) next
  auc <- get_max_roc_auc(seurat_obj, pair[1], pair[2])
  if (auc < auto_merge_threshold) {
    seurat_obj <- SetIdent(seurat_obj,
                           cells = WhichCells(seurat_obj, idents = pair[2]),
                           value = pair[1])
    message("Auto‑merged ", pair[2], " into ", pair[1], " (AUC = ", round(auc, 3), ")")
  }
}

seurat_obj$seurat_clusters <- Idents(seurat_obj)
message("Final number of clusters: ", length(levels(seurat_obj)))

# Custom Plot for Figure 3D
library(ggplot2)
library(dplyr)

# Extract t-SNE coordinates and metadata
tsne_coords <- Embeddings(seurat_obj, reduction = "tsne")
plot_df <- data.frame(
  tSNE_1 = tsne_coords[, 1],
  tSNE_2 = tsne_coords[, 2],
  cluster = Idents(seurat_obj),                
  sample  = seurat_obj$stage,                  
  region  = seurat_obj$region                  
)

plot_df$cluster <- factor(plot_df$cluster,
                          levels = levels(plot_df$cluster),
                          labels = paste0("C", as.numeric(levels(plot_df$cluster)) + 1))

plot_df <- plot_df %>%
  mutate(
    shape_id = case_when(
      !is.na(region) & region == "r1" ~ "r1_53d",
      !is.na(region) & region == "r2" ~ "r2_53d",
      !is.na(region) & region == "r3" ~ "r3_58d",
      !is.na(region) & region == "r4" ~ "r4_58d",
      sample == "33 days" ~ "33d",
      sample == "35 days" ~ "35d",
      sample == "37 days" ~ "37d",
      sample == "41 days" ~ "41d",
      sample == "65 days" ~ "65d",
      TRUE ~ NA_character_
    ),
    is_whole = !grepl("^r[1-4]", shape_id)
  )

shape_mapping <- c(
  "r1_53d" = 22,   # square
  "r2_53d" = 23,   # diamond
  "r3_58d" = 21,   # circle
  "r4_58d" = 24,   # triangle point up
  "33d"    = 21,   # circle
  "35d"    = 24,   # triangle
  "37d"    = 25,   # upside‑down triangle
  "41d"    = 23,   # diamond
  "65d"    = 22    # square
)

n_clusters <- length(levels(plot_df$cluster))
cluster_colors <- scales::hue_pal()(n_clusters)

df_whole <- plot_df %>% filter(is_whole)
df_micro <- plot_df %>% filter(!is_whole)

# Create the plot
p_custom <- ggplot() +
  geom_point(data = df_whole,
             aes(x = tSNE_1, y = tSNE_2, shape = shape_id,
                 fill = cluster, color = cluster),
             size = 2.5, stroke = 0.8) +
  geom_point(data = df_micro,
             aes(x = tSNE_1, y = tSNE_2, shape = shape_id, color = cluster),
             fill = "white", size = 2.5, stroke = 0.8) +
  
  scale_shape_manual(values = shape_mapping) +
  scale_fill_manual(values = cluster_colors) +   
  scale_color_manual(values = cluster_colors) + 
  
  theme_minimal() +
  theme(
    legend.position = "right"
  ) +
  guides(
    shape = guide_legend(title = "Sample", override.aes = list(color = "black")),
    color = guide_legend(title = "Cluster"),
    fill  = guide_legend(title = "Cluster") 
  ) +
  ggtitle("Organoid single‑cell clusters (Fig. 3D)")

print(p_custom)

p_custom <- recordPlot()

saveRDS(p_custom, file = "Fig3D_plot.rds")


# Export cluster assignments
cluster_assignments <- Idents(seurat_obj)

levels(cluster_assignments) <- paste0("C", as.numeric(levels(cluster_assignments)) + 1)

cells_by_cluster <- split(names(cluster_assignments), cluster_assignments)

max_len <- max(sapply(cells_by_cluster, length))

cluster_df <- as.data.frame(lapply(cells_by_cluster, function(x) {
  c(x, rep(NA, max_len - length(x)))
}))

write.csv(cluster_df, file = "Organoidclusters.csv", row.names = FALSE, na = "")

message("Cluster assignments written to Organoidclusters.csv")


# Gene expression feature plots on t‑SNE (marker genes)

marker_genes <- c("FOXG1", "OTX2", "RSPO2", "DCN", "ASPM", "LIN28A", "MYT1L", "NEUROD6")

feature_plots <- lapply(marker_genes, function(gene) {
  if (gene %in% rownames(seurat_obj)) {
    FeaturePlot(seurat_obj,
                features = gene,
                reduction = "tsne",
                cols = c("grey90", "red"),
                order = TRUE) +
      ggtitle(gene) +
      theme_void() +                               # removes axes, grid, background
      theme(plot.title = element_text(face = "italic", hjust = 0.5),
            legend.position = "none")              # no individual legends
  } else NULL
})

feature_plots <- Filter(Negate(is.null), feature_plots)

# Create dummy plot to extract legend
dummy <- FeaturePlot(seurat_obj,
                     features = marker_genes[1],
                     reduction = "tsne",
                     cols = c("grey90", "red"),
                     order = TRUE) +
  scale_color_gradient(low = "grey90", high = "red",
                       breaks = range(FetchData(seurat_obj, vars = marker_genes[1])[,1]),
                       labels = c("Low", "High")) +
  guides(color = guide_colorbar(title = NULL)) +
  theme_void() +                                   # keep legend background clean
  theme(legend.position = "right")

library(cowplot)
shared_legend <- get_legend(dummy)

# Combine plots and shared legend
combined_plots <- wrap_plots(feature_plots, ncol = 4)
final_plot <- plot_grid(combined_plots, shared_legend,
                        ncol = 2, rel_widths = c(0.9, 0.1))

print(final_plot)