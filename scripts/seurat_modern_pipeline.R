library(Seurat)
library(dplyr)
library(ggplot2)
library(plotly)
library(htmlwidgets)

set.seed(1234)

# Main Seurat reanalysis script

counts_raw <- read.table(
  "GSE75140_hOrg.fetal.master.data.frame.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

# Matrix comes in as cells x genes, so flip it for Seurat
counts <- t(counts_raw)

seurat_obj <- CreateSeuratObject(
  counts = counts,
  min.cells = 3,
  min.features = 200
)

# Basic Filtering 
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 500 &
           nFeature_RNA < 7000
)

# No extra mitochondrial filtering here since I was mainly keeping the same broad workflow 
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)

# PCA + clustering
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# 0.5 was the most reasonable middle ground, lower merged too much, higher started splitting things too hard
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)

png("outputs/UMAP_2D_clusters.png", width = 1000, height = 800)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
dev.off()

# Also saved a 3D version just to check cluster separation another way
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

# Marker genes by cluster
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

# Simple condition labels from the sample names
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
