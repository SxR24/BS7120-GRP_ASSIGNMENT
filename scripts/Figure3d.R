#Figure 3D
# Load the libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Load the data
raw <- read.table(
  "/home/ngk6/Desktop/BS7120 Steered research project (Group)/GSE75140_hOrg.fetal.master.data.frame.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = ""
)

#Split metadata vs expression
META_COLS <- c(
  "cell_id", "Cell", "CellID", "barcode",
  "sample", "Sample", "donor", "Donor",
  "day", "Day", "stage", "Stage", "timepoint", "Timepoint",
  "batch", "Batch", "replicate", "Replicate",
  "species", "Species"
)

meta_cols <- intersect(colnames(raw), META_COLS)
gene_cols <- setdiff(colnames(raw), meta_cols)

if (length(gene_cols) < 500) {
  stop("Too few gene columns detected. Likely you need to fix META_COLS based on your file.")
}

meta <- raw[, meta_cols, drop = FALSE]
expr <- raw[, gene_cols, drop = FALSE]

# Converting expression to numeric matrix
expr_num <- as.matrix(sapply(expr, function(x) suppressWarnings(as.numeric(x))))

#Replace NA with 0 
expr_num[is.na(expr_num)] <- 0

# Give cells rownames
cell_ids <- NULL
if ("cell_id" %in% colnames(meta)) cell_ids <- meta$cell_id
if (is.null(cell_ids) && "CellID" %in% colnames(meta)) cell_ids <- meta$CellID
if (is.null(cell_ids) && "Cell" %in% colnames(meta)) cell_ids <- meta$Cell
if (is.null(cell_ids)) cell_ids <- paste0("Cell_", seq_len(nrow(expr_num)))

rownames(expr_num) <- make.unique(as.character(cell_ids))
colnames(expr_num) <- gene_cols

# Seurat Analysis
counts <- t(expr_num)

#Creating the Seurat Object
seurat_obj <- CreateSeuratObject(
  counts = counts,
  meta.data = meta,
  min.cells = 3,
  min.features = 200
)

# QC Metrics 
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

#Cleaning 
nf <- seurat_obj$nFeature_RNA
nc <- seurat_obj$nCount_RNA
pm <- seurat_obj$percent.mt

nf_lo <- unname(quantile(nf, 0.01))
nf_hi <- unname(quantile(nf, 0.99))
nc_lo <- unname(quantile(nc, 0.01))
nc_hi <- unname(quantile(nc, 0.99))

#Removal of mitochondria 
pm_hi <- 10

seurat_obj <- subset(
  seurat_obj,
  subset =
    nFeature_RNA >= nf_lo & nFeature_RNA <= nf_hi &
    nCount_RNA   >= nc_lo & nCount_RNA   <= nc_hi &
    percent.mt   <= pm_hi
)

#Remove genes expressed in extremely few cells (extra cleanup)
# seurat_obj <- seurat_obj[Matrix::rowSums(seurat_obj@assays$RNA@counts > 0) >= 3, ]

#Normalization + HVGs
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

#Scaling
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("nCount_RNA", "percent.mt"))

#Dim Reduction
seurat_obj <- RunPCA(seurat_obj, npcs = 30, features = VariableFeatures(seurat_obj), verbose = FALSE)

pcs_to_use <- 1:20
seurat_obj <- FindNeighbors(seurat_obj, dims = pcs_to_use, verbose = FALSE)

#Clustering
target_k <- 11
res_grid <- seq(0.2, 2.0, by = 0.02)

best_res <- NA
best_diff <- Inf

for (r in res_grid) {
  tmp <- FindClusters(seurat_obj, resolution = r, verbose = FALSE)
  k <- length(unique(Idents(tmp)))
  d <- abs(k - target_k)
  
  if (d < best_diff) {
    best_diff <- d
    best_res <- r
  }
  if (k == target_k) {
    best_res <- r
    break
  }
}

seurat_obj <- FindClusters(seurat_obj, resolution = best_res, verbose = FALSE)
cat("Chosen resolution:", best_res, " | clusters:", length(unique(Idents(seurat_obj))), "\n")

#Rename clusters 
lvl <- levels(Idents(seurat_obj))
new_ids <- paste0("c", seq_along(lvl))
names(new_ids) <- lvl
seurat_obj <- RenameIdents(seurat_obj, new_ids)

# t-SNE plot
set.seed(42)
seurat_obj <- RunTSNE(seurat_obj, dims = pcs_to_use, perplexity = 30, verbose = FALSE)

# Plotting
tp_candidates <- c("day", "Day", "stage", "Stage", "timepoint", "Timepoint", "sample", "Sample")
tp_col <- intersect(tp_candidates, colnames(seurat_obj@meta.data))
tp_col <- if (length(tp_col) > 0) tp_col[1] else NULL

p <- DimPlot(
  seurat_obj,
  reduction = "tsne",
  label = TRUE,
  pt.size = 0.8,
  shape.by = tp_col
) + ggtitle("t-SNE with 11 clusters (c1–c11)")

print(p)
ggsave("tsne_11clusters_clean.png", plot = p, width = 7.5, height = 6.5, dpi = 300)

# MARKERS
markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(top10, "top10_markers_11clusters.csv", row.names = FALSE)
