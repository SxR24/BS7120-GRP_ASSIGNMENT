# ------------------------------------------------------------
# Modern Seurat Reanalysis of Camp et al. (2015)
# Dataset: GSE75140
# ------------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(plotly)
library(htmlwidgets)

set.seed(1234)

# ------------------------------------------------------------
# 1. Load Data
# ------------------------------------------------------------

counts_raw <- read.table(
  "GSE75140_hOrg.fetal.master.data.frame.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

# Transpose (dataset is cells x genes; Seurat needs genes x cells)
counts <- t(counts_raw)

seurat_obj <- CreateSeuratObject(
  counts = counts,
  min.cells = 3,
  min.features = 200
)

# ------------------------------------------------------------
# 2. Basic Filtering
# ------------------------------------------------------------

seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 500 &
           nFeature_RNA < 7000
)

# ------------------------------------------------------------
# 3. SCTransform Normalisation
# ------------------------------------------------------------

seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)

# ------------------------------------------------------------
# 4. PCA
# ------------------------------------------------------------

seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

# ------------------------------------------------------------
# 5. Clustering (Resolution 0.5)
# ------------------------------------------------------------

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# ------------------------------------------------------------
# 6. 2D UMAP
# ------------------------------------------------------------

seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)

png("outputs/UMAP_2D_clusters.png", width = 1000, height = 800)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
dev.off()

# ------------------------------------------------------------
# 7. 3D UMAP
# ------------------------------------------------------------

seurat_obj <- RunUMAP(seurat_obj, dims = 1:15, n.components = 3)
umap_3d <- Embeddings(seurat_obj, "umap")

plot3D <- plot_ly(
  x = umap_3d[,1],
  y = umap_3d[,2],
  z = umap_3d[,3],
  type = "scatter3d",
  mode = "markers",
  color = as.factor(Idents(seurat_obj)),
  marker = list(size = 3)
)

htmlwidgets::saveWidget(
  plot3D,
  "outputs/Final_UMAP_3D_resolution0.5.html"
)

# ------------------------------------------------------------
# 8. Marker Identification
# ------------------------------------------------------------

markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(
  top_markers,
  "outputs/Top10_markers_per_cluster.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Condition Comparison (Fetal vs Organoid)
# ------------------------------------------------------------

seurat_obj$condition <- ifelse(
  grepl("hOrg", colnames(seurat_obj)),
  "Organoid",
  "Fetal"
)

cluster_condition_table <- table(
  Idents(seurat_obj),
  seurat_obj$condition
)

write.csv(
  cluster_condition_table,
  "outputs/Cluster_vs_Condition_table.csv"
)

# ------------------------------------------------------------
# End of Script
# ------------------------------------------------------------
