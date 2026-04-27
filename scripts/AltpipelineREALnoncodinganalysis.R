# =============================================================================
# Alt-pipeline noncoding QC steps
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(Matrix)
library(leiden)
library(biomaRt)
library(scran)
library(SingleCellExperiment)

set.seed(42)

# Inputs

counts_file <- "D:/GroupA_University2026_Project/ALT_pipeline_data/counts.txt"

counts <- tryCatch({
  read.table(counts_file, header = TRUE, row.names = 1, 
             sep = "\t", check.names = FALSE)
}, error = function(e) {
  read.table(counts_file, header = TRUE, row.names = 1, 
             sep = ",", check.names = FALSE)
})

cat("=== COUNT MATRIX INFO ===\n")
cat("Dimensions (genes x cells):", dim(counts), "\n")
cat("First 5 row names (genes):", head(rownames(counts), 5), "\n")

# Convert to sparse matrix
counts <- as.matrix(counts)
counts <- as(counts, "dgCMatrix")

# Detect ensemble genes and map to symbol IDs
is_ensembl <- all(grepl("^ENS", rownames(counts)[1:min(10, nrow(counts))]))

if (is_ensembl) {
  cat("\nDetected Ensembl IDs. Mapping to gene symbols using biomaRt...\n")
  
  # Human Ensembl to Gene Symbol mapping
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  gene_map <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = rownames(counts),
    mart = ensembl
  )
  
  gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
  rownames(gene_map) <- gene_map$ensembl_gene_id
  
  original_ids <- rownames(counts)
  new_names <- gene_map[original_ids, "external_gene_name"]
  new_names[is.na(new_names) | new_names == ""] <- original_ids[is.na(new_names) | new_names == ""]
  
  dup_symbols <- duplicated(new_names) | duplicated(new_names, fromLast = TRUE)
  new_names[dup_symbols] <- paste0(new_names[dup_symbols], "_", original_ids[dup_symbols])
  
  rownames(counts) <- new_names
  cat("Mapped", sum(!is.na(gene_map$external_gene_name)), "genes to symbols.\n")
} else {
  cat("\nRow names are not Ensembl IDs. Proceeding with original names.\n")
}

# Create seurat obejct

seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "ALT_pipeline",
  min.cells = 10,
  min.features = 200
)

cat("\n=== SEURAT OBJECT INITIAL ===\n")
print(seurat_obj)
cat("Cells before QC:", ncol(seurat_obj), "\n")

# Mitochondrial gene detection

mito_pattern <- "^MT-"
test_mito <- grep(mito_pattern, rownames(seurat_obj), value = TRUE, ignore.case = TRUE)

if (length(test_mito) == 0) {
  cat("\nNo 'MT-' genes found. Using Ensembl mitochondrial IDs...\n")
  
  # Human mitochondrial Ensembl IDs (GRCh38)
  # These are the 13 protein-coding, 2 rRNA, and 22 tRNA genes
  mito_ensembl_ids <- c(
    "ENSG00000198888", "ENSG00000198763", "ENSG00000198804",
    "ENSG00000212907", "ENSG00000198786", "ENSG00000198695",
    "ENSG00000198727", "ENSG00000198840", "ENSG00000210082",
    "ENSG00000198938", "ENSG00000198899", "ENSG00000198712",
    "ENSG00000210100", "ENSG00000210107", "ENSG00000210112",
    "ENSG00000210117", "ENSG00000210127", "ENSG00000210135",
    "ENSG00000210140", "ENSG00000210144", "ENSG00000210151",
    "ENSG00000210154", "ENSG00000210156", "ENSG00000210164",
    "ENSG00000210174", "ENSG00000210176", "ENSG00000210184",
    "ENSG00000210191", "ENSG00000210194", "ENSG00000210195",
    "ENSG00000210196", "ENSG00000209082", "ENSG00000211459"
  )
  
  # Find which of these are in the dataset
  mito_genes_present <- intersect(mito_ensembl_ids, rownames(seurat_obj))
  
  if (length(mito_genes_present) > 0) {
    # Create a custom pattern that matches these exact IDs
    mito_pattern <- paste0("^(", paste(mito_genes_present, collapse = "|"), ")$")
    cat("Found", length(mito_genes_present), "mitochondrial genes in dataset.\n")
  } else {
    warning("No mitochondrial genes found. Percent.mt will be zero.")
    mito_pattern <- "^NONEXISTENT$"
  }
}

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mito_pattern)

# QC Visualization

# Generate individual violin plots
vln_list <- VlnPlot(seurat_obj,
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    combine = FALSE, pt.size = 0.1)

vln_list[[1]] <- vln_list[[1]] + 
  labs(x = "", y = "Genes Detected", title = "Genes per Cell")
vln_list[[2]] <- vln_list[[2]] + 
  labs(x = "", y = "Total Reads", title = "Reads per Cell")
vln_list[[3]] <- vln_list[[3]] + 
  labs(x = "", y = "% Mitochondrial", title = "Mitochondrial %")

p1 <- wrap_plots(vln_list, ncol = 3)

p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  labs(x = "Total Reads", y = "% Mitochondrial", title = "Mitochondrial Content")

p3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  labs(x = "Total Reads", y = "Genes Detected", title = "Library Complexity")

combined_plot <- p1 / (p2 | p3) +
  plot_annotation(title = "Quality Control (Before Filtering)",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

# Save and display
ggsave("QC_before_filtering.png", plot = combined_plot, width = 14, height = 10)
print(combined_plot)

cat("\n=== QC STATISTICS ===\n")
cat("nFeature_RNA quantiles:\n")
print(quantile(seurat_obj$nFeature_RNA, probs = c(0.01, 0.05, 0.5, 0.95, 0.99)))
cat("nCount_RNA quantiles:\n")
print(quantile(seurat_obj$nCount_RNA, probs = c(0.01, 0.05, 0.5, 0.95, 0.99)))
cat("percent.mt quantiles:\n")
print(quantile(seurat_obj$percent.mt, probs = c(0.01, 0.05, 0.5, 0.95, 0.99)))

# QC thresholds

nFeature_lower <- 2000     # Fluidigm C1 typically captures >2000 genes per cell
nFeature_upper <- 12000    # Based on 99th percentile
nCount_lower <- 50000      # Read counts
nCount_upper <- 1200000    # Keep large cells, remove only extreme outlineres
percent_mt_max <- 10       # Mitochondrial filter

seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA > nFeature_lower & 
                       nFeature_RNA < nFeature_upper &
                       nCount_RNA > nCount_lower & 
                       nCount_RNA < nCount_upper &
                       percent.mt < percent_mt_max)

cat("\nCells remaining after QC:", ncol(seurat_obj), "\n")

# Normalization

# Convert Seurat to SingleCellExperiment
sce <- as.SingleCellExperiment(seurat_obj)

sce <- computeSumFactors(sce, clusters = quickCluster(sce))

sf <- sizeFactors(sce)

# Apply normalization to Seurat object
raw_counts <- GetAssayData(seurat_obj, layer = "counts")

# Compute normalized data: log1p( counts / size_factor )
norm_data <- log1p(t(t(raw_counts) / sf))

seurat_obj[["RNA"]]$data <- norm_data

seurat_obj$size_factor <- sf

print(head(seurat_obj[["RNA"]]$data[1:5, 1:3]))

# After PCA
saveRDS(seurat_obj, file = "seurat_preprocessed.rds")
scaled_matrix <- GetAssayData(seurat_obj, layer = "scale.data")
saveRDS(scaled_matrix, file = "D:/GroupA_University2026_Project/OUTPUTS/seurat_matrix.rds")
