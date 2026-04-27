# =====================================================================================
# Alt noncoding pipeline: subpopulation assignment to fetal and organoid clusters
# (Per-cell correlations, NPC subanalysis, and t‑SNE for organoid)
# =====================================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Matrix)
library(SingleCellExperiment)
library(scran)
library(leiden)
library(biomaRt)
library(pheatmap)
library(future)
library(stringr)
library(cowplot)


plan("multisession", workers = 4, gc = TRUE)
options(future.globals.maxSize = 2000 * 1024^2)
set.seed(12345, kind = "L'Ecuyer-CMRG")   # reproducible parallel RNG

# Load pre‑processed Seurat object and remove interneurons samples
seurat_obj <- readRDS("seurat_preprocessed.rds")
seurat_obj <- subset(seurat_obj, cells = setdiff(colnames(seurat_obj), 
                                                 c("SRR2967628", "SRR2967645", "SRR2967696", "SRR2967704", "SRR2967725")))

# Load metadata
meta <- read.csv("metadata.csv", stringsAsFactors = FALSE)
region_meta <- read.csv("Organoidlabelmetadata.csv", stringsAsFactors = FALSE)

# Add tissue type metadata (Fetal vs Organoid)
metadata_file <- "metadata.csv"
meta <- read.csv(metadata_file, stringsAsFactors = FALSE)

cell_names <- colnames(seurat_obj)
run_ids <- meta$Run

if (all(cell_names %in% run_ids)) {
  sample_type <- meta$tissue[match(cell_names, meta$Run)]
}

cell_barcodes <- colnames(organoid_obj)
meta_match <- meta[match(cell_barcodes, meta$cell), "stage"]
region_match <- region_meta[match(cell_barcodes, region_meta$cell), "region"]

organoid_obj$sample <- meta_match
organoid_obj$region <- region_match

seurat_obj$tissue_type <- ifelse(sample_type == "Fetal neocortex", "Fetal", "Organoid")
cat("Sample type distribution:\n")
print(table(seurat_obj$tissue_type))

# Split into separate Fetal and Organoid objects
fetal_obj <- subset(seurat_obj, subset = tissue_type == "Fetal")
organoid_obj <- subset(seurat_obj, subset = tissue_type == "Organoid")
cat("Fetal cells:", ncol(fetal_obj), "\n")
cat("Organoid cells:", ncol(organoid_obj), "\n")

# Define a standardised clustering workflow
#    Uses variable genes for initial PCA, selects genes with high loadings,
#    re‑runs PCA on the selected set, then performs graph‑based clustering
#    and UMAP visualisation.
analyze_subset <- function(obj, name, 
                           n_pcs_for_selection = 15,
                           loading_threshold = 0.2,
                           genes_per_pc = 50,
                           n_pcs_for_clustering = 5) {
  cat("\n--- Analyzing", name, "---\n")
  if (ncol(obj) < 50) {
    warning(sprintf("Only %d cells in %s. Skipping clustering.", ncol(obj), name))
    return(obj)
  }
  
  DefaultAssay(obj) <- "RNA"
  
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
  var_genes <- VariableFeatures(obj)
  cat("Using", length(var_genes), "variable genes for initial PCA.\n")
  
  obj <- ScaleData(obj, features = var_genes, verbose = FALSE)
  obj <- RunPCA(obj, features = var_genes, npcs = max(n_pcs_for_selection, n_pcs_for_clustering), verbose = FALSE)
  
  loadings <- Loadings(obj, reduction = "pca")
  selected_genes <- c()
  
  for (i in 1:n_pcs_for_selection) {
    pc_loadings <- loadings[, i]
    strong_genes <- names(pc_loadings)[abs(pc_loadings) > loading_threshold]
    if (length(strong_genes) == 0) next
    strong_genes <- strong_genes[order(abs(pc_loadings[strong_genes]), decreasing = TRUE)]
    top_genes <- head(strong_genes, genes_per_pc)
    selected_genes <- union(selected_genes, top_genes)
    cat("PC", i, ":", length(strong_genes), "genes above threshold, added", length(top_genes), "\n")
  }
  
  cat("Total selected genes from initial PCA:", length(selected_genes), "\n")
  if (length(selected_genes) < 50) warning("Fewer than 50 genes selected; consider adjusting thresholds.")
  
  # Re‑run PCA on the refined gene set
  VariableFeatures(obj) <- selected_genes
  obj <- ScaleData(obj, features = selected_genes, verbose = FALSE)
  obj <- RunPCA(obj, features = selected_genes, npcs = n_pcs_for_clustering, verbose = FALSE)
  
  # Elbow plot for diagnostics
  ElbowPlot(obj, ndims = n_pcs_for_clustering)
  ggsave(paste0("Elbow_", name, "_loading_selected.png"), width = 8, height = 6)
  print(ElbowPlot(obj, ndims = n_pcs_for_clustering))
  
  # Graph‑based clustering and UMAP
  n_pcs_use <- min(n_pcs_for_clustering, ncol(obj) - 1)
  
  obj <- FindNeighbors(obj, dims = 1:n_pcs_use, k.param = 20, prune.SNN = 1/15)
  obj <- FindClusters(obj, resolution = 2, algorithm = 4)
  
  obj <- RunUMAP(obj, dims = 1:n_pcs_use)
}

DimPlot(obj, reduction = "umap")
# Apply the clustering workflow to Fetal and Organoid objects
fetal_obj <- analyze_subset(fetal_obj, "Fetal")
organoid_obj <- analyze_subset(organoid_obj, "Organoid")

# t‑SNE visualisation for organoid
tsne_success <- FALSE
organoid_obj <- tryCatch({
  obj <- RunTSNE(organoid_obj, 
                 dims = 1:min(10, ncol(organoid_obj) - 1),
                 perplexity = min(30, floor((ncol(organoid_obj) - 1) / 3 - 1)),
                 check_duplicates = FALSE)
  tsne_success <- TRUE
  obj
}, error = function(e) {
  warning("t‑SNE failed: ", e$message, ". No Fig. 3D plot will be generated. 
            Adjust perplexity or number of cells and re‑run.")
  organoid_obj
})

if (tsne_success) {
  # Add stage and region metadata (ensure metadata.csv has 'stage' and 'region' columns)
  organoid_obj$stage  <- meta$stage[match(colnames(organoid_obj), meta$Run)]
  organoid_obj$region <- meta$region[match(colnames(organoid_obj), meta$Run)]
  
  # Extract t‑SNE coordinates
  coords <- Embeddings(organoid_obj, reduction = "tsne")
  plot_df <- data.frame(
    tSNE_1 = coords[, 1],
    tSNE_2 = coords[, 2],
    cluster = Idents(organoid_obj),
    stage   = organoid_obj$stage,
    region  = organoid_obj$region,
    stringsAsFactors = FALSE
  )
  
  plot_df$cluster <- factor(plot_df$cluster,
                            levels = levels(plot_df$cluster),
                            labels = paste0("C", as.numeric(levels(plot_df$cluster)) + 1))
  
  # Create shape identifiers
  plot_df <- plot_df %>%
    mutate(
      shape_id = case_when(
        !is.na(region) & region == "r1" ~ "r1_53d",
        !is.na(region) & region == "r2" ~ "r2_53d",
        !is.na(region) & region == "r3" ~ "r3_58d",
        !is.na(region) & region == "r4" ~ "r4_58d",
        stage == "33 days" ~ "33d",
        stage == "35 days" ~ "35d",
        stage == "37 days" ~ "37d",
        stage == "41 days" ~ "41d",
        stage == "65 days" ~ "65d",
        TRUE ~ NA_character_
      ),
      is_whole = !grepl("^r[1-4]", shape_id)
    )
  
  # Shape mapping (as in the original Fig. 3D)
  shape_mapping <- c(
    "r1_53d" = 22,   
    "r2_53d" = 23,   
    "r3_58d" = 21,   
    "r4_58d" = 24,   
    "33d"    = 21,   
    "35d"    = 24,
    "37d"    = 25,   
    "41d"    = 23,   
    "65d"    = 22    
  )
  
  # Cluster colours
  n_clusters <- length(levels(plot_df$cluster))
  cluster_colors <- scales::hue_pal()(n_clusters)
  
  # Split data for different point aesthetics
  df_whole <- plot_df %>% filter(is_whole)
  df_micro <- plot_df %>% filter(!is_whole)
  
  # Build the plot
  p_fig3d <- ggplot() +
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
    labs(x = "tSNE 1", y = "tSNE 2",
         title = "Organoid single‑cell clusters (Fig. 3D)") +
    theme_minimal() +
    theme(legend.position = "right") +
    guides(shape = guide_legend(title = "Sample", override.aes = list(color = "black")),
           color = guide_legend(title = "Cluster"),
           fill  = guide_legend(title = "Cluster"))
  
  print(p_fig3d)
  ggsave("Fig3D_organoid_tsne.png", p_fig3d, width = 10, height = 7, dpi = 300)
  saveRDS(p_fig3d, file = "Fig3D_plot.rds")
  
} else {
  message("Fig. 3D plot skipped because t‑SNE could not be computed.")
}

# Correlate clusters with bulk microdissected zones (W13)
bulk_file <- "Old_GSE38805_human_FPKM_matrix_log2_symbols.txt"
bulk_matrix <- read.table(bulk_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

w13_cols <- grep("^w13_", colnames(bulk_matrix), value = TRUE)
bulk_zones <- bulk_matrix[, w13_cols, drop = FALSE]
zone_names <- sub("^w13_", "", w13_cols)
colnames(bulk_zones) <- zone_names
bulk_zones <- as.matrix(bulk_zones)

common_genes <- intersect(rownames(fetal_obj), rownames(bulk_zones))
if (length(common_genes) < 500) warning("Few common genes between single-cell and bulk data.")

# Function: assign each cluster to the best‑correlated zone
assign_zones <- function(obj, bulk_mat, common_genes, cor_threshold = 0.2) {
  cluster_avg <- AggregateExpression(obj, assays = "RNA", slot = "data",
                                     features = common_genes,
                                     group.by = "seurat_clusters",
                                     return.seurat = FALSE)$RNA
  cluster_avg <- as.matrix(cluster_avg)
  bulk_sub <- bulk_mat[rownames(cluster_avg), , drop = FALSE]
  
  cor_mat <- cor(cluster_avg, bulk_sub, method = "spearman")
  
  max_cor <- apply(cor_mat, 1, max)
  best_zone <- colnames(cor_mat)[apply(cor_mat, 1, which.max)]
  names(best_zone) <- rownames(cor_mat)
  
  weak_clusters <- names(max_cor)[max_cor < cor_threshold]
  if (length(weak_clusters) > 0) {
    warning("Clusters with max correlation < ", cor_threshold, ": ", 
            paste(weak_clusters, collapse = ", "))
    best_zone[weak_clusters] <- "Unassigned"
  }
  
  obj$predicted_zone <- best_zone[as.character(obj$seurat_clusters)]
  return(list(obj = obj, cor_mat = cor_mat))
}

# Run zone assignment on both objects
fetal_res <- assign_zones(fetal_obj, bulk_zones, common_genes)
fetal_obj <- fetal_res$obj
cor_fetal <- fetal_res$cor_mat

organoid_res <- assign_zones(organoid_obj, bulk_zones, common_genes)
organoid_obj <- organoid_res$obj
cor_organoid <- organoid_res$cor_mat

# Per‑cell correlation with microdissected zones
cat("Computing per-cell correlations with cortical zones...\n")
expr_mat <- LayerData(fetal_obj, assay = "RNA", layer = "data")
common_genes <- intersect(rownames(expr_mat), rownames(bulk_zones))
if (length(common_genes) < 50) warning("Fewer than 50 common genes; correlations may be unreliable.")

expr_sub <- as.matrix(expr_mat[common_genes, , drop = FALSE])
bulk_sub <- as.matrix(bulk_zones[common_genes, , drop = FALSE])

cell_cor <- cor(expr_sub, bulk_sub, method = "spearman")
cell_cor_scaled <- t(scale(t(cell_cor)))

fetal_obj$CP_cor_scaled   <- cell_cor_scaled[, "CP"]
fetal_obj$ISVZ_cor_scaled <- cell_cor_scaled[, "ISVZ"]
fetal_obj$OSVZ_cor_scaled <- cell_cor_scaled[, "OSVZ"]
fetal_obj$VZ_cor_scaled   <- cell_cor_scaled[, "VZ"]

# Combined UMAP of zone correlations
zone_plot_list <- lapply(c("CP", "ISVZ", "OSVZ", "VZ"), function(zn) {
  FeaturePlot(fetal_obj, features = paste0(zn, "_cor_scaled"), reduction = "umap") +
    ggtitle(paste("Fetal", zn, "correlation")) +
    theme(legend.position = "none")
})
combined_zones <- wrap_plots(zone_plot_list, ncol = 2)
print(combined_zones)

# Heatmap (only if ≤5000 cells)
if (ncol(fetal_obj) < 5000) {
  pheatmap::pheatmap(cell_cor_scaled,
                     show_rownames = FALSE,
                     treeheight_row = 0,
                     treeheight_col = 0,
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     main = "Correlation with microdissected zones")
}

# =============================================================================
# NPC SUBPOPULATION DISCOVERY (small‑n clustering on NPC cells)
# =============================================================================
npc_file <- "noncoding_Sample_lists_Neurons_vs_NPCs.tsv"
if (file.exists(npc_file)) {
  cat("Loading NPC/Neuron metadata...\n")
  npc_meta <- read.table(npc_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  npc_samples <- na.omit(npc_meta$NPCs)
  npc_samples <- npc_samples[nchar(npc_samples) > 0]
  npc_cells <- intersect(npc_samples, colnames(fetal_obj))
  cat("Identified", length(npc_cells), "NPC cells.\n")
  
  if (length(npc_cells) >= 10) {
    npc_obj <- subset(fetal_obj, cells = npc_cells)
    
    DefaultAssay(npc_obj) <- "RNA"
    npc_obj <- NormalizeData(npc_obj, verbose = FALSE)
    npc_obj <- FindVariableFeatures(npc_obj, selection.method = "vst", nfeatures = 1000, verbose = FALSE)
    npc_obj <- ScaleData(npc_obj, verbose = FALSE)
    
    n_pcs <- min(4, ncol(npc_obj) - 1)
    npc_obj <- RunPCA(npc_obj, features = VariableFeatures(npc_obj), npcs = n_pcs, verbose = FALSE)
    print(ElbowPlot(npc_obj, ndims = n_pcs) + ggtitle("NPC subset - PCA Elbow"))
    
    k_val <- min(5, ncol(npc_obj) - 1)
    npc_obj <- FindNeighbors(npc_obj, dims = 1:n_pcs, k.param = k_val, verbose = FALSE)
    npc_obj <- FindClusters(npc_obj, resolution = 2, algorithm = 4, verbose = FALSE)
    
    n_neighbors_umap <- min(15, ncol(npc_obj) - 1)
    npc_obj <- RunUMAP(npc_obj, dims = 1:n_pcs, n.neighbors = n_neighbors_umap, verbose = FALSE)
    
    p_npc_clusters <- DimPlot(npc_obj, reduction = "umap", label = TRUE, label.size = 5) +
      ggtitle(paste0("NPC Subclusters (n = ", ncol(npc_obj), " cells)")) +
      theme_minimal()
    print(p_npc_clusters)
    ggsave("NPC_subclusters_UMAP.png", p_npc_clusters, width = 8, height = 6, dpi = 300)
    saveRDS(npc_obj, file = "npc_subsetevil_annotated.rds")
    
    if (length(unique(npc_obj$seurat_clusters)) > 1) {
      npc_markers <- FindAllMarkers(npc_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
      write.csv(npc_markers, "NPC_subcluster_markers.csv", row.names = FALSE)
      top_markers <- npc_markers %>% group_by(cluster) %>% slice_head(n = 10)
      print(top_markers)
      
      top_genes <- unique(top_markers$gene)
      if (length(top_genes) > 1) {
        p_heat <- DoHeatmap(npc_obj, features = top_genes, size = 3) + 
          ggtitle("Top NPC Subcluster Markers")
        ggsave("NPC_marker_heatmap.png", p_heat, width = 10, height = 8, dpi = 300)
      }
    }
    cat("NPC subclustering complete. Found", length(unique(npc_obj$seurat_clusters)), "clusters.\n")
  } else {
    warning("Fewer than 10 NPC cells identified; clustering not meaningful.")
  }
}

# Subcluster cells within the VZ zone only
subcluster_cells_safe <- function(obj, cell_names, name_prefix) {
  if (length(cell_names) < 10) {
    cat("Skipping", name_prefix, "- too few cells (", length(cell_names), ")\n")
    return(NULL)
  }
  sub_obj <- subset(obj, cells = cell_names)
  n_cells <- ncol(sub_obj)
  max_pcs <- min(10, n_cells - 1, nrow(sub_obj) - 1)
  cat("Subclustering", name_prefix, "with", n_cells, "cells. Using", max_pcs, "PCs.\n")
  
  sub_obj <- FindVariableFeatures(sub_obj, nfeatures = 1000)
  sub_obj <- ScaleData(sub_obj)
  sub_obj <- RunPCA(sub_obj, npcs = max_pcs, verbose = FALSE)
  sub_obj <- FindNeighbors(sub_obj, dims = 1:max_pcs)
  sub_obj <- FindClusters(sub_obj, resolution = 0.8, algorithm = 4)
  sub_obj <- RunUMAP(sub_obj, dims = 1:max_pcs)
  
  print(DimPlot(sub_obj, reduction = "umap", label = TRUE) +
          ggtitle(paste(name_prefix, "Subclusters (n =", n_cells, ")")))
  return(sub_obj)
}

vz_cells <- colnames(fetal_obj)[fetal_obj$predicted_zone == "VZ"]
vz_cells <- vz_cells[!is.na(vz_cells)]

if (length(vz_cells) > 0) {
  vz_sub <- subcluster_cells_safe(fetal_obj, cell_names = vz_cells, name_prefix = "Fetal_VZ")
  if (!is.null(vz_sub)) {
    fetal_obj$VZ_subcluster <- NA
    fetal_obj$VZ_subcluster[colnames(vz_sub)] <- as.character(vz_sub$seurat_clusters)
  }
}

# Save the final annotated objects
saveRDS(fetal_obj, file = "fetal_noncodingannotated.rds")
saveRDS(organoid_obj, file = "organoid_noncodingannotated.rds")

# =============================================================================
# Feature plots of marker genes on organoid t‑SNE / UMAP
# =============================================================================
marker_genes <- c("FOXG1", "OTX2", "RSPO2", "DCN", "ASPM", "LIN28A", "MYT1L", "NEUROD6")

if (tsne_success && "tsne" %in% Reductions(organoid_obj)) {
  plot_obj <- organoid_obj
  reduction_use <- "tsne"
} else {
  plot_obj <- organoid_obj
  reduction_use <- ifelse("umap" %in% Reductions(plot_obj), "umap", "pca")
}
cat("Using", reduction_use, "for organoid marker feature plots.\n")

feature_plots <- lapply(marker_genes, function(gene) {
  if (gene %in% rownames(plot_obj)) {
    FeaturePlot(plot_obj, features = gene, reduction = reduction_use,
                cols = c("grey90", "red"), order = TRUE) +
      ggtitle(gene) + theme_void() +
      theme(plot.title = element_text(face = "italic", hjust = 0.5),
            legend.position = "none")
  } else {
    warning("Gene ", gene, " not found in dataset.")
    NULL
  }
})
feature_plots <- Filter(Negate(is.null), feature_plots)

if (length(feature_plots) > 0) {
  dummy <- FeaturePlot(plot_obj, features = marker_genes[1], reduction = reduction_use,
                       cols = c("grey90", "red"), order = TRUE) +
    scale_color_gradient(low = "grey90", high = "red",
                         breaks = range(FetchData(plot_obj, vars = marker_genes[1])[,1]),
                         labels = c("Low", "High")) +
    guides(color = guide_colorbar(title = NULL)) +
    theme_void() + theme(legend.position = "right")
  shared_legend <- get_legend(dummy)
  final_plot <- plot_grid(wrap_plots(feature_plots, ncol = 4), shared_legend,
                          ncol = 2, rel_widths = c(0.9, 0.1))
  print(final_plot)
} else {
  cat("No marker genes found. Skipping feature plots.\n")
}

# Fetal marker feature plots (UMAP)
fetal_markers <- c("PAX6", "VIM", "CCNB1", "MKI67", "NEUROD6", "NEUROD4", "MYT1L", "MEF2C")
reduction_fetal <- ifelse("umap" %in% Reductions(fetal_obj), "umap", "pca")

fetal_feature_plots <- lapply(fetal_markers, function(gene) {
  if (gene %in% rownames(fetal_obj)) {
    FeaturePlot(fetal_obj, features = gene, reduction = reduction_fetal,
                cols = c("grey90", "red"), order = TRUE) +
      ggtitle(gene) + theme_void() +
      theme(plot.title = element_text(face = "italic", hjust = 0.5),
            legend.position = "none")
  } else NULL
})
fetal_feature_plots <- Filter(Negate(is.null), fetal_feature_plots)

if (length(fetal_feature_plots) > 0) {
  first_marker <- fetal_markers[fetal_markers %in% rownames(fetal_obj)][1]
  dummy_fetal <- FeaturePlot(fetal_obj, features = first_marker, reduction = reduction_fetal,
                             cols = c("grey90", "red"), order = TRUE) +
    scale_color_gradient(low = "grey90", high = "red",
                         breaks = range(FetchData(fetal_obj, vars = first_marker)[,1]),
                         labels = c("Low", "High")) +
    guides(color = guide_colorbar(title = NULL)) +
    theme_void() + theme(legend.position = "right")
  shared_legend_fetal <- get_legend(dummy_fetal)
  final_plot_fetal <- plot_grid(wrap_plots(fetal_feature_plots, ncol = 4), shared_legend_fetal,
                                ncol = 2, rel_widths = c(0.9, 0.1))
  print(final_plot_fetal)
}

# NPC subset marker feature plots
if (exists("npc_obj") && ncol(npc_obj) > 10) {
  npc_markers <- c("PAX6", "VIM", "CCNB1", "MKI67", "MYT1L", "MEF2C", "EOMES", "ASPM")
  reduction_npc <- ifelse("umap" %in% Reductions(npc_obj), "umap", "pca")
  
  npc_feature_plots <- lapply(npc_markers, function(gene) {
    if (gene %in% rownames(npc_obj)) {
      FeaturePlot(npc_obj, features = gene, reduction = reduction_npc,
                  cols = c("grey90", "red"), order = TRUE) +
        ggtitle(gene) + theme_void() +
        theme(plot.title = element_text(face = "italic", hjust = 0.5),
              legend.position = "none")
    } else NULL
  })
  npc_feature_plots <- Filter(Negate(is.null), npc_feature_plots)
  
  if (length(npc_feature_plots) > 0) {
    first_marker_npc <- npc_markers[npc_markers %in% rownames(npc_obj)][1]
    dummy_npc <- FeaturePlot(npc_obj, features = first_marker_npc, reduction = reduction_npc,
                             cols = c("grey90", "red"), order = TRUE) +
      scale_color_gradient(low = "grey90", high = "red",
                           breaks = range(FetchData(npc_obj, vars = first_marker_npc)[,1]),
                           labels = c("Low", "High")) +
      guides(color = guide_colorbar(title = NULL)) +
      theme_void() + theme(legend.position = "right")
    shared_legend_npc <- get_legend(dummy_npc)
    final_plot_npc <- plot_grid(wrap_plots(npc_feature_plots, ncol = 4), shared_legend_npc,
                                ncol = 2, rel_widths = c(0.9, 0.1))
    print(final_plot_npc)
  }
}

# =============================================================================
# Export PC gene lists from fetal PCA
# =============================================================================
if (!"pca" %in% Reductions(fetal_obj)) stop("PCA reduction not found.")
all_loadings <- Loadings(fetal_obj, reduction = "pca")

# Define which PCs and signs to extract
my_pc_selections <- list(
  list(PC = 1, sign = +1, label = "PC1 cor"),
  list(PC = 1, sign = -1, label = "PC1 anti"),
  list(PC = 2, sign = +1, label = "PC2 cor"),
  list(PC = 2, sign = -1, label = "PC2 anti"),
  list(PC = 3, sign = -1, label = "PC3 anti"),
  list(PC = 4, sign = -1, label = "PC4 anti"),
  list(PC = 1, sign = +1, label = "PC1 rg cor"),
  list(PC = 2, sign = +1, label = "PC2 rg cor"),
  list(PC = 2, sign = -1, label = "PC2 rg anti")
)

get_top_genes <- function(loadings, pc_num, sign, n_genes = 100) {
  if (pc_num > ncol(loadings)) return(character(0))
  pc_loadings <- loadings[, pc_num]
  idx <- if (sign == 1) pc_loadings > 0 else pc_loadings < 0
  sorted_genes <- names(sort(abs(pc_loadings[idx]), decreasing = TRUE))
  head(sorted_genes, n_genes)
}

pc_lists <- lapply(my_pc_selections, function(sel) {
  get_top_genes(loadings = all_loadings, pc_num = sel$PC, sign = sel$sign, n_genes = 100)
})
names(pc_lists) <- sapply(my_pc_selections, `[[`, "label")

# Pad and save as CSV
max_len <- max(sapply(pc_lists, length))
pc_lists_padded <- lapply(pc_lists, function(x) {
  if (length(x) < max_len) c(x, rep(NA, max_len - length(x))) else x
})
write.csv(as.data.frame(pc_lists_padded), 
          "D:/GroupA_University2026_Project/AltnoncodingPCA.csv",
          row.names = FALSE, na = "")
message("PC gene lists exported.")

# =============================================================================
# Per‑cell correlation with FACS‑purified aRG, bRG, N (GSE65000)
# =============================================================================
bulk_file_aRG <- "Old_GSE65000_hsa_fpkm_matrix_log2.txt"
if (file.exists(bulk_file_aRG)) {
  cat("Loading GSE65000 for aRG/bRG/N correlations...\n")
  bulk_mat_arg <- read.table(bulk_file_aRG, header = TRUE, row.names = 1,
                             sep = "\t", check.names = FALSE)
  
  aRG_cols <- grep("^Hsa_aRG_", colnames(bulk_mat_arg), value = TRUE)
  bRG_cols <- grep("^Hsa_bRG_", colnames(bulk_mat_arg), value = TRUE)
  N_cols   <- grep("^Hsa_N_",   colnames(bulk_mat_arg), value = TRUE)
  
  bulk_arg_avg <- data.frame(
    aRG = rowMeans(bulk_mat_arg[, aRG_cols, drop = FALSE], na.rm = TRUE),
    bRG = rowMeans(bulk_mat_arg[, bRG_cols, drop = FALSE], na.rm = TRUE),
    N   = rowMeans(bulk_mat_arg[, N_cols,   drop = FALSE], na.rm = TRUE)
  )
  rownames(bulk_arg_avg) <- rownames(bulk_mat_arg)
  
  expr_mat <- LayerData(fetal_obj, assay = "RNA", layer = "data")
  common_genes_arg <- intersect(rownames(expr_mat), rownames(bulk_arg_avg))
  if (length(common_genes_arg) >= 50) {
    expr_sub_arg <- as.matrix(expr_mat[common_genes_arg, , drop = FALSE])
    bulk_sub_arg <- as.matrix(bulk_arg_avg[common_genes_arg, c("aRG", "bRG", "N")])
    
    cell_cor_arg <- cor(expr_sub_arg, bulk_sub_arg, method = "spearman")
    cell_cor_arg_scaled <- t(scale(t(cell_cor_arg)))
    colnames(cell_cor_arg_scaled) <- c("aRG", "bRG", "N")
    
    fetal_obj$aRG_cor_scaled <- cell_cor_arg_scaled[, "aRG"]
    fetal_obj$bRG_cor_scaled <- cell_cor_arg_scaled[, "bRG"]
    fetal_obj$N_cor_scaled   <- cell_cor_arg_scaled[, "N"]
    
    if (ncol(fetal_obj) < 5000) {
      pheatmap::pheatmap(cell_cor_arg_scaled,
                         show_rownames = FALSE,
                         treeheight_row = 0, treeheight_col = 0,
                         color = colorRampPalette(c("blue", "white", "red"))(100),
                         main = "Correlation with FACS-purified cell types")
    }
  } else {
    warning("Fewer than 50 common genes for aRG/bRG/N.")
  }
} else {
  warning("GSE65000 file not found.")
}

# =============================================================================
# Per‑cell zone assignment based on max scaled correlation
# =============================================================================
zone_cor_cols <- c("CP_cor_scaled", "ISVZ_cor_scaled", "OSVZ_cor_scaled", "VZ_cor_scaled")
zone_cor_mat <- fetal_obj@meta.data[, zone_cor_cols, drop = FALSE]
max_zone_idx <- apply(zone_cor_mat, 1, which.max)
zone_names_clean <- c("CP", "iSVZ", "oSVZ", "VZ")
assigned_zone_per_cell <- zone_names_clean[max_zone_idx]
assigned_zone_per_cell[apply(zone_cor_mat, 1, function(x) all(is.na(x)))] <- NA
fetal_obj$predicted_zone <- assigned_zone_per_cell

cat("Per‑cell zone distribution:\n")
print(table(fetal_obj$predicted_zone, useNA = "always"))


# Export cell assignments (zone + cell type) to text file

cat("Creating altpipeline_noncoding_cell_assignments file...\n")
zone_out <- fetal_obj$predicted_zone

# Cell type from aRG/bRG/N correlations (if available)
if (all(c("aRG_cor_scaled", "bRG_cor_scaled", "N_cor_scaled") %in% colnames(fetal_obj@meta.data))) {
  type_cor_mat <- fetal_obj@meta.data[, c("aRG_cor_scaled", "bRG_cor_scaled", "N_cor_scaled")]
  max_type_idx <- apply(type_cor_mat, 1, which.max)
  type_map <- c("aRG", "bRG", "neuron")
  assigned_type <- type_map[max_type_idx]
  assigned_type[apply(type_cor_mat, 1, function(x) all(is.na(x)))] <- NA
} else {
  assigned_type <- rep(NA, ncol(fetal_obj))
}

write.table(data.frame(Sample = colnames(fetal_obj),
                       Assigned_Zone = zone_out,
                       Assigned_CellType = assigned_type,
                       stringsAsFactors = FALSE),
            file = "altpipeline_noncoding_cell_assignments.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat("Assignments exported.\n")

# =============================================================================
# Custom UMAP with NPC subcluster overlays
# =============================================================================
group_cols <- c("AP2" = "#95F297", "AP1" = "#2C8334",
                "BP1" = "#15B5B5", "BP2" = "#9DFFFF",
                "N1"  = "#8EC9F5", "N2"  = "#4589EF", "N3"  = "#0545BB")

fetal_obj$plot_group <- as.character(fetal_obj$seurat_clusters)

# Assign broad neuron and progenitor groups
fetal_obj$plot_group[fetal_obj$seurat_clusters %in% c(5, 3)] <- "N3"
fetal_obj$plot_group[fetal_obj$seurat_clusters == 2]      <- "N1"
fetal_obj$plot_group[fetal_obj$seurat_clusters == 6]      <- "N2"

# Overlay detailed NPC subcluster identities
if (exists("npc_obj") && ncol(npc_obj) >= 10) {
  npc_cells <- colnames(npc_obj)
  subclust <- as.character(npc_obj$seurat_clusters)
  names(subclust) <- npc_cells
  
  sub_to_type <- c(
    "5" = "BP1", "4" = "AP1", "6" = "AP1",
    "3" = "AP2", "1" = "AP2", "2" = "BP2"
  )
  for (cell in npc_cells) {
    sub_id <- subclust[cell]
    if (sub_id %in% names(sub_to_type)) fetal_obj$plot_group[cell] <- sub_to_type[sub_id]
  }
  cat("NPC overlay applied.\n")
}

# Build colour vector and plot
all_groups <- unique(fetal_obj$plot_group)
full_cols <- group_cols[names(group_cols) %in% all_groups]
remaining_groups <- setdiff(all_groups, names(full_cols))
if (length(remaining_groups) > 0) {
  hues <- scales::hue_pal()(length(remaining_groups))
  names(hues) <- remaining_groups
  full_cols <- c(full_cols, hues)
}

umap_coords <- Embeddings(fetal_obj, reduction = "umap")
plot_data <- data.frame(UMAP_1 = umap_coords[, 1],
                        UMAP_2 = umap_coords[, 2],
                        group = fetal_obj$plot_group)

p_manual <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, fill = group)) +
  geom_point(shape = 21, colour = "black", stroke = 0.2, size = 2.4) +
  scale_fill_manual(values = full_cols) +
  ggtitle("Fetal Neocortex Subpopulations") +
  theme_minimal() +
  theme(legend.position = "right", axis.title = element_text(size = 12))
print(p_manual)

# Quick check for cluster 4 markers
cluster4_markers <- FindMarkers(fetal_obj, ident.1 = 4, min.pct = 0.25)
head(cluster4_markers, 20)

# Export final subpopulation assignments
subpop_df <- data.frame(sample_id = colnames(fetal_obj),
                        group = fetal_obj$plot_group,
                        stringsAsFactors = FALSE)
write.csv(subpop_df, "Altnoncodingsubpopulations.csv", row.names = FALSE, quote = FALSE)
cat("Subpopulation assignments exported. Group counts:\n")
print(table(subpop_df$group))