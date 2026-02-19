# ===============================
# Seurat RNA-seq Reanalysis
# BS7120 Group Assignment
# ===============================

library(Seurat)
library(dplyr)
library(ggplot2)
library(plotly)
library(htmlwidgets)

# -------------------------------
# 1. Load Data
# -------------------------------

counts <- read.table("GSE75140_hOrg.fetal.master.data.frame.txt",
                     header = TRUE,
                     row.names = 1)

# -------------------------------
# 2. Create Seurat Object
# -------------------------------

seurat_obj <- CreateSeuratObject(
  counts = counts,
  min.cells = 3,
  min.features = 200
)

# -------------------------------
# 3. QC Metrics
# -------------------------------

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^MT-"
)

seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 500 & nFeature_RNA < 7000
)

# -------------------------------
# 4. Normalization (SCTransform)
# -------------------------------

seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)

# -------------------------------
# 5. PCA
# -------------------------------

seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

# -------------------------------
# 6. Neighbors + Clustering
# -------------------------------

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# -------------------------------
# 7. UMAP 2D
# -------------------------------

seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)

p2d <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)

ggsave("UMAP_2D_clusters.png",
       plot = p2d,
       width = 7,
       height = 6)

# -------------------------------
# 8. UMAP 3D
# -------------------------------

seurat_obj <- RunUMAP(seurat_obj,
                      dims = 1:15,
                      n.components = 3)

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
  "Final_UMAP_3D_resolution0.5.html"
)

# -------------------------------
# 9. Marker Detection
# -------------------------------

markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(top_markers,
          "Top10_markers_per_cluster.csv",
          row.names = FALSE)

# -------------------------------
# 10. Condition Annotation
# -------------------------------

seurat_obj$condition <- ifelse(
  grepl("hOrg", colnames(seurat_obj)),
  "Organoid",
  "Fetal"
)

cluster_condition_table <- table(
  Idents(seurat_obj),
  seurat_obj$condition
)

write.csv(cluster_condition_table,
          "Cluster_vs_Condition_table.csv")
