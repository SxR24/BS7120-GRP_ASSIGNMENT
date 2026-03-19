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
library(Matrix)
library(patchwork)
library(ggplot2)

# Monocle pipeline

# Import cell metadata and expression matrix; no gene metadata as not required

cell_metadata<- read_delim("cell_metadata_1.tsv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)
expression_matrix <- read_delim("expression_matrix2.tsv", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE)
gene_metadata <- read_delim("gene_metadata.tsv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

# Adjust dfs and convert expression_matrix into a matrix format

expression_matrix <- expression_matrix %>% tibble::column_to_rownames(., var = names(.)[1])
cell_metadata <- cell_metadata %>% tibble::column_to_rownames(., var = names(.)[1])
gene_metadata <- gene_metadata %>% tibble::column_to_rownames(., var = names(.)[1])
mat <- as.matrix(expression_matrix)


# Monocles main class; cell_data_set

cds <- new_cell_data_set(mat,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata) 

# Cells are to be arranged based off of genes used to classify cells from Fig1C; here are the cells

genes_from_fig1c <- c("PAX6", "GLI3", "VIM","SOX2","MYT1L","TBR1","NEUROD6","DCX",
                      "FAT3","CDH6","CDH7","SNCA","ROBO2",
                      "DMD","MKI67","NRCAM", "NRP1", "POU3F2", "ELMO1", "SEMA5A", "PDPN", "KIF11", "NEUROD4", "EOMES",
                      "HES1","HES6","ASPM","BCL11B"
                 ) 

# Genes selected as of Fig2B

fig2b <- c("SOX2","EOMES","MYT1L")

# Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, use_genes = genes_from_fig1c)
# Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds)
# Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)
# Step 4: Cluster the cells
cds <- cluster_cells(cds)
# Step 5: Learn a graph
cds <- learn_graph(cds)
# Step 6: Order cells (Starts at AP)
cds <- order_cells(cds)


a <- plot_cells(cds, color_cells_by = "CellType",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           cell_size = 2)

b <- plot_cells(cds, color_cells_by = "Zone",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           cell_size = 2)
c <- plot_cells(cds, genes = "SOX2",
                label_cell_groups=FALSE,
                label_leaves=TRUE,
                label_branch_points=TRUE,
                graph_label_size=1.5,
                cell_size = 2) 
d <- plot_cells(cds, genes = "EOMES",
                label_cell_groups=FALSE,
                label_leaves=TRUE,
                label_branch_points=TRUE,
                graph_label_size=1.5,
                cell_size = 2) + theme(legend.position = "none")
e <- plot_cells(cds, genes = "MYT1L",
                label_cell_groups=FALSE,
                label_leaves=TRUE,
                label_branch_points=TRUE,
                graph_label_size=1.5,
                cell_size = 2) + theme(legend.position = "none")
a / (b+c+d+e)

save_monocle_objects(cds=cds, directory_path='monocle_3_cds', comment='Version1 - 18/03')

## Monocle Objects - Please read for LAMP server!

# Best way to load up the 'monocle_3_cds' is via:
# cds <- load_monocle_objects(directory_path='monocle_3_cds')
#
# May have to change file path, but doing so will initialise the cds in the current session,
# mitigating the need for the raw files and processing on the server.
# From here, just make graphs as shown above (via plot_cells())
# For reference, valid 'colouring' would be:
# 'Zone', 'Cluster', 'CellType', 'stage'
# 
# It may be useful to have a drop-down for genes; not entirely sure how shiny works but
# see above for a list of genes used, it may be nice to visualize different gene expressions

