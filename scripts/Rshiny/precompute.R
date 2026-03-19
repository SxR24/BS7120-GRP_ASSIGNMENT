# precompute_gse75140.R
# ------------------------------------------------------------
# Run once to preprocess the dataset and save fast-loading RDS files
# ------------------------------------------------------------

library(data.table)
library(Seurat)
library(Matrix)
library(dplyr)

# -----------------------------
# Settings
# -----------------------------
DATA_FILE <- "C:\\Users\\Sohil Ananth\\Downloads\\t3_testing\\GSE75140_hOrg.fetal.master.data.frame.txt"
SEED_USE  <- 42

PCS_TO_USE <- 1:20
TARGET_K   <- 11
RES_GRID   <- seq(0.2, 2.0, by = 0.02)

QC_Q_LO <- 0.01
QC_Q_HI <- 0.99
PM_HI   <- 10

# -----------------------------
# Helpers
# -----------------------------
read_gse75140 <- function(path) {
  df <- fread(
    file = path,
    sep = "\t",
    header = TRUE,
    data.table = FALSE,
    quote = "",
    fill = TRUE
  )
  
  colnames(df) <- gsub('"', "", colnames(df))
  
  if ("cell_id" %in% colnames(df)) {
    df$cell_id <- gsub('"', "", df$cell_id)
  }
  
  if ("species" %in% colnames(df)) {
    df$species <- gsub('"', "", df$species)
  }
  
  if (!all(c("cell_id", "species") %in% colnames(df))) {
    stop("Expected columns 'cell_id' and 'species' not found.")
  }
  
  df
}

split_meta_expr <- function(raw) {
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
    stop("Too few gene columns detected. Check input format.")
  }
  
  list(
    meta = raw[, meta_cols, drop = FALSE],
    gene_cols = gene_cols
  )
}

choose_resolution_for_k <- function(seu_obj, res_grid, target_k) {
  best_res  <- NA_real_
  best_diff <- Inf
  
  for (r in res_grid) {
    tmp <- FindClusters(seu_obj, resolution = r, verbose = FALSE)
    k <- length(unique(Idents(tmp)))
    d <- abs(k - target_k)
    
    if (d < best_diff) {
      best_diff <- d
      best_res  <- r
    }
    
    if (k == target_k) {
      best_res <- r
      break
    }
  }
  
  list(best_res = best_res, best_diff = best_diff)
}

rename_clusters_c <- function(seu_obj) {
  lvl <- levels(Idents(seu_obj))
  new_ids <- paste0("c", seq_along(lvl))
  names(new_ids) <- lvl
  RenameIdents(seu_obj, new_ids)
}

get_embed_df <- function(seu_obj, reduction = c("tsne", "umap", "pca")) {
  reduction <- match.arg(reduction)
  
  emb <- Embeddings(seu_obj, reduction = reduction)
  df <- as.data.frame(emb)
  df$cell_id <- rownames(df)
  
  md <- seu_obj@meta.data[, "species", drop = FALSE]
  md$cell_id <- rownames(md)
  
  df <- dplyr::left_join(df, md, by = "cell_id")
  df$cluster <- as.character(Idents(seu_obj)[df$cell_id])
  
  df
}

# -----------------------------
# Main pipeline
# -----------------------------
set.seed(SEED_USE)

cat("Reading dataset...\n")
raw <- read_gse75140(DATA_FILE)

cat("Splitting metadata and expression...\n")
sp <- split_meta_expr(raw)
meta <- sp$meta
gene_cols <- sp$gene_cols

expr <- raw[, gene_cols, drop = FALSE]

cat("Converting expression to numeric matrix...\n")
expr_num <- as.matrix(sapply(expr, function(x) suppressWarnings(as.numeric(x))))
expr_num[is.na(expr_num)] <- 0

cell_ids <- NULL
if ("cell_id" %in% colnames(meta)) cell_ids <- meta$cell_id
if (is.null(cell_ids) && "CellID" %in% colnames(meta)) cell_ids <- meta$CellID
if (is.null(cell_ids) && "Cell" %in% colnames(meta)) cell_ids <- meta$Cell
if (is.null(cell_ids)) cell_ids <- paste0("Cell_", seq_len(nrow(expr_num)))

rownames(expr_num) <- make.unique(as.character(cell_ids))
colnames(expr_num) <- gene_cols

cat("Building sparse count matrix...\n")
counts <- Matrix(t(expr_num), sparse = TRUE)

cat("Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(
  counts = counts,
  meta.data = meta,
  min.cells = 3,
  min.features = 200
)

cat("Calculating mitochondrial percentage...\n")
mt_genes <- grep("^MT-", rownames(seurat_obj), value = TRUE)

if (length(mt_genes) > 0) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
} else {
  seurat_obj$percent.mt <- 0
}

cat("Applying QC filters...\n")
nf <- seurat_obj$nFeature_RNA
nc <- seurat_obj$nCount_RNA

nf_lo <- unname(quantile(nf, QC_Q_LO))
nf_hi <- unname(quantile(nf, QC_Q_HI))
nc_lo <- unname(quantile(nc, QC_Q_LO))
nc_hi <- unname(quantile(nc, QC_Q_HI))

seurat_obj <- subset(
  seurat_obj,
  subset =
    nFeature_RNA >= nf_lo & nFeature_RNA <= nf_hi &
    nCount_RNA   >= nc_lo & nCount_RNA   <= nc_hi &
    percent.mt   <= PM_HI
)

cat("Normalizing data...\n")
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

cat("Finding variable features...\n")
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

cat("Scaling data...\n")
seurat_obj <- ScaleData(
  seurat_obj,
  vars.to.regress = c("nCount_RNA", "percent.mt"),
  verbose = FALSE
)

cat("Running PCA...\n")
seurat_obj <- RunPCA(
  seurat_obj,
  npcs = 30,
  features = VariableFeatures(seurat_obj),
  verbose = FALSE
)

cat("Finding neighbors...\n")
seurat_obj <- FindNeighbors(
  seurat_obj,
  dims = PCS_TO_USE,
  verbose = FALSE
)

cat("Choosing clustering resolution...\n")
res_pick <- choose_resolution_for_k(seurat_obj, RES_GRID, TARGET_K)

cat("Best resolution:", res_pick$best_res, "\n")
cat("Best diff:", res_pick$best_diff, "\n")

cat("Clustering cells...\n")
seurat_obj <- FindClusters(
  seurat_obj,
  resolution = res_pick$best_res,
  verbose = FALSE
)

cat("Renaming clusters to c1..cN...\n")
seurat_obj <- rename_clusters_c(seurat_obj)

cat("Running tSNE...\n")
set.seed(SEED_USE)
seurat_obj <- RunTSNE(
  seurat_obj,
  dims = PCS_TO_USE,
  perplexity = min(30, floor((ncol(seurat_obj) - 1) / 3)),
  verbose = FALSE
)

cat("Running UMAP...\n")
seurat_obj <- RunUMAP(
  seurat_obj,
  dims = PCS_TO_USE,
  verbose = FALSE
)

cat("Extracting embeddings...\n")
emb_tsne <- get_embed_df(seurat_obj, "tsne")
emb_umap <- get_embed_df(seurat_obj, "umap")
emb_pca  <- get_embed_df(seurat_obj, "pca")

cat("Saving RDS files...\n")
saveRDS(seurat_obj, "gse75140_seurat_ready.rds")
saveRDS(emb_tsne, "emb_tsne.rds")
saveRDS(emb_umap, "emb_umap.rds")
saveRDS(emb_pca, "emb_pca.rds")

cat("Done.\n")
cat("Saved:\n")
cat("- gse75140_seurat_ready.rds\n")
cat("- emb_tsne.rds\n")
cat("- emb_umap.rds\n")
cat("- emb_pca.rds\n")