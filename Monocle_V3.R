## MONOCLE SETUP

# For installing monocle 3, ensure the latest version of BiocManager has been installed
# https://cole-trapnell-lab.github.io/monocle3/docs/installation/
# See code below; taken from above site.

remove.packages("monocle") # Ignore if no previous experience with monocle
remove.packages("BiocManager") # Ignore if no previous experience with BiocManager

# Note, if there are issues removing BiocManager *and* the below code doesn't resolve as intended
# navigate to /lib/R/site-library and remove BiocManager via console: "sudo rm -rf BiocManager"


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.21") # Previous version 1.30.27

# Rmpfr fails; Header file mpfr.h not found... 
# Configuration failed for 'units'; libudunits2.so not found
# Configuration failed for 'gdtools; cairo-ft.h not found
# Configuration failed for 'sf'; gdal-config not found or not exectuable
# Configuration failed for 'exactextractr'; geos-config not found or not executable
# Loading failed for 'stars'; found sf but wrong version (>= 1.0.19 required)
# Potential Other Errors - TL;DR lots of errors but so far no problems

# Caution - Above install takes ~ 2.5 hours.

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'ggrastr')) # More dependencies

library(devtools)


remotes::install_github(c("bnprks/BPCells/r")) # Required by Monocle 3

# If the above doesn't work, there is a good chance it is due to not finding 'hdf5.h'.
# I recommend adding the following lines to your .bashrc file and retrying: 
# export CPATH=/opt/homebrew/include
# export LIBRARY_PATH=/opt/homebrew/lib
# Thanks to: https://github.com/bnprks/BPCells/issues/6#issuecomment-1700435556


devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")


## Start of Monocle3 pipeline

library(monocle3)
library(readr)
library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratWrappers)
library(Matrix) # Not used for Seurat - only if using 'raw' data


# Seurat Pipeline
PATH <- "/INSERT/PATH/HERE"
pbmc.data <- Read10X(data.dir = PATH)
# Initialize the Seurat object with the raw (non-normalized data).

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Monocle from Seurat Object {Recommended}
# Please see https://stuartlab.org/signac/articles/monocle for more detail

pbmc.cds <- as.cell_data_set(pbmc)
pbmc.cds <- cluster_cells(pbmc.cds, reduction_method = "UMAP")
pbmc.cds <- learn_graph(pbmc.cds)
pbmc.cds <- order_cells(pbmc.cds)

plot_cells(
  cds = pbmc.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)


# Monocle from raw data {Potential Issues} 
# Works best from pre-existing Seurat Object

matrix <- Matrix::readMM("PATH/matrix.mtx")

barcodes <- read_csv("PATH/barcodes.tsv", 
                     col_names = FALSE)
genes <- read_table("PATH/genes.tsv", 
                    col_names = FALSE)
# May not be necessary, check the file first.
colnames(genes) <- c("GENCODE", "gene_short_name")

cds <- new_cell_data_set(matrix,
                         cell_metadata = barcodes,
                         gene_metadata = genes) # Monocles main class; cell_data_set


# Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100) # Note - Takes a moment (~1 minute)
# Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds, alignment_group = "batch")
# Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

# Step 4: Cluster the cells
cds <- cluster_cells(cds)
# Step 5: Learn a graph
cds <- learn_graph(cds)

# Step 6: Order cells
cds <- order_cells(cds)

plot_cells(cds)