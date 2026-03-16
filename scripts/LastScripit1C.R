# Load the required libraries 
library(data.table) 
library(matrixStats) 
library(pheatmap)
# # Define the path to the gene expression matrix ( must be changed)
expr_file <- "/home/nhsas1/Desktop/2026_UoL_prj_gene_matrix_f.tsv"
# # Load the expression matrix
expr <- fread(expr_file)
# Extract gene names and remove the gene column from the matrix
gene_names <- expr$gene
expr$gene <- NULL
# Convert the data table into a numeric matrix
expr_mat <- as.matrix(expr)

#ensure numeric

expr_mat <- apply(expr_mat, 2, as.numeric)
# Assign gene names as row names of the matrix
rownames(expr_mat) <- gene_names

# inspect matrix
dim(expr_mat)
expr_mat[1:5,1:5]
summary(expr_mat)
# Load metadata containing SRR IDs, sample type, stage, and chip information 
meta_table <- read.csv("/home/nhsas1/Desktop/Rtest/SRR_PaperID_Metadata.csv")
# Identify SRR IDs corresponding to fetal cells present in the matrix
fetal_srr <- intersect(
  meta_table$SRR[meta_table$type == "fetal"],
  colnames(expr_mat)
)
# Check number of fetal cells recovered
length(fetal_srr)
# Subset the expression matrix to retain only fetal cells
expr_fetal <- expr_mat[, fetal_srr, drop = FALSE]
# Confirm dimensions of fetal-only
dim(expr_fetal)
# Verify that the endothelial marker gene PECAM1 exists
if(!"PECAM1" %in% rownames(expr_fetal)) {
  stop("PECAM1 not found in expression matrix")
}
# Transpose matrix so rows represent cells and columns represent genes
expr_cells <- t(expr_fetal)

#inspect PECAM1 distribution

sort(expr_cells[, "PECAM1"], decreasing = TRUE)[1:10]
# Identify the cell with highest PECAM1 expression (endothelial cell)
endo_cell <- rownames(expr_cells)[which.max(expr_cells[, "PECAM1"])]
# Confirm PECAM1 expression in the identified endothelial cell
endo_cell
expr_cells[endo_cell, "PECAM1"]

#remove endothelial cell

expr_no_endo <- expr_fetal[, colnames(expr_fetal) != endo_cell]
# Confirm new matrix dimensions after endothelial removal
dim(expr_no_endo)
# Transpose matrix so rows represent cells for marker scoring
expr_cells <- t(expr_no_endo)
# Define known interneuron marker genes
interneuron_markers <- c(
  "DLX1","DLX2","DLX5","DLX6",
  "GAD1","GAD2","LHX6"
)

#keep only markers present in dataset

interneuron_markers <- interneuron_markers[
  interneuron_markers %in% colnames(expr_cells)
]
# Check which marker genes are available
interneuron_markers
length(interneuron_markers)

#calculate interneuron score

interneuron_score <- rowSums(
  expr_cells[, interneuron_markers, drop = FALSE] > 1
)
# Inspect distribution of marker scores across cells
table(interneuron_score)
# Identify cells expressing multiple interneuron markers
interneuron_cells <- names(
  interneuron_score[interneuron_score >= 2]
)
# Check number of detected interneuron cells
length(interneuron_cells)

#inspect marker expression

expr_cells[interneuron_cells, interneuron_markers]

#remove interneurons

expr_clean <- expr_no_endo[
  , !(colnames(expr_no_endo) %in% interneuron_cells)
]
# Confirm dimensions after removing endothelial and interneurons
dim(expr_clean)
# Count number of cells expressing each gene (>1 FPKM)
expr_count <- rowSums(expr_clean > 1)
# Inspect distribution of gene detection counts
summary(expr_count)
# Keep genes expressed in more than two cells
expr_expr <- expr_clean[expr_count > 2, ]
# Inspect matrix dimensions after expression filtering
dim(expr_expr)
# Compute variance of each gene across cells
gene_var <- rowVars(expr_expr)
# Inspect distribution of gene variances
summary(gene_var)
# Keep highly variable genes (variance > 0.5)
expr_var <- expr_expr[gene_var > 0.5, ]
# Inspect matrix after variance filtering
dim(expr_var)
# prepare PCA input
# Transpose matrix for PCA (cells × genes)
pca_input <- t(expr_var)

# remove genes with NA
pca_input <- pca_input[, colSums(is.na(pca_input)) == 0]

# verify matrix
sum(is.na(pca_input))
sum(!is.finite(pca_input))

# run PCA
pca_res <- prcomp(
  pca_input,
  center = TRUE,
  scale. = FALSE
)
# Inspect variance explained by principal components
summary(pca_res)

# Prepare PC score matrix
#
# Extract cell scores from PCA. These represent each cell's
# position along each principal component.

pc_scores <- as.data.frame(pca_res$x)
# Inspect PCA score matrix dimensions and first components
dim(pc_scores)
head(pc_scores[, 1:3])

# Remove genes with missing names before correlation analysis
sum(is.na(rownames(expr_var)))
expr_var <- expr_var[!is.na(rownames(expr_var)), ]
sum(is.na(rownames(expr_var)))

#Compute gene-PC correlations
#   For each gene, correlate its expression across cells with
#   the cell PC scores. This gives interpretable gene-PC
#   associations that can be thresholded at 0.2.

gene_pc_cor <- data.frame(
  gene = rownames(expr_var),
  PC1  = apply(expr_var, 1, function(g) cor(as.numeric(g), pc_scores$PC1)),
  PC2  = apply(expr_var, 1, function(g) cor(as.numeric(g), pc_scores$PC2)),
  PC3  = apply(expr_var, 1, function(g) cor(as.numeric(g), pc_scores$PC3))
)
# Preview gene-PC correlation results
head(gene_pc_cor)

# Extract PC1–PC3 correlated and anticorrelated genes
#   absolute loading/correlation > 0.2
#   maximum 50 genes per PC direction

# # Define correlation threshold and maximum genes per PC direction
loading_threshold <- 0.2
max_genes <- 50

# PC1 correlated genes

pc1_cor_df <- gene_pc_cor[gene_pc_cor$PC1 > loading_threshold, ]
pc1_cor_df <- pc1_cor_df[order(pc1_cor_df$PC1, decreasing = TRUE), ]
pc1_cor <- head(pc1_cor_df$gene, max_genes)

# PC1 anticorrelated genes
pc1_anti_df <- gene_pc_cor[gene_pc_cor$PC1 < -loading_threshold, ]
pc1_anti_df <- pc1_anti_df[order(pc1_anti_df$PC1, decreasing = FALSE), ]
pc1_anti <- head(pc1_anti_df$gene, max_genes)

# PC2 correlated genes

pc2_cor_df <- gene_pc_cor[gene_pc_cor$PC2 > loading_threshold, ]
pc2_cor_df <- pc2_cor_df[order(pc2_cor_df$PC2, decreasing = TRUE), ]
pc2_cor <- head(pc2_cor_df$gene, max_genes)


# PC2 anticorrelated genes

pc2_anti_df <- gene_pc_cor[gene_pc_cor$PC2 < -loading_threshold, ]
pc2_anti_df <- pc2_anti_df[order(pc2_anti_df$PC2, decreasing = FALSE), ]
pc2_anti <- head(pc2_anti_df$gene, max_genes)

# PC3 correlated genes


pc3_cor_df <- gene_pc_cor[gene_pc_cor$PC3 > loading_threshold, ]
pc3_cor_df <- pc3_cor_df[order(pc3_cor_df$PC3, decreasing = TRUE), ]
pc3_cor <- head(pc3_cor_df$gene, max_genes)


# PC3 anticorrelated genes

pc3_anti_df <- gene_pc_cor[gene_pc_cor$PC3 < -loading_threshold, ]
pc3_anti_df <- pc3_anti_df[order(pc3_anti_df$PC3, decreasing = FALSE), ]
pc3_anti <- head(pc3_anti_df$gene, max_genes)


# Check number of genes selected for each PC direction
length(pc1_cor)
length(pc1_anti)
length(pc2_cor)
length(pc2_anti)
length(pc3_cor)
length(pc3_anti)

# inspect example genes
head(pc1_cor)
head(pc1_anti)
head(pc2_cor)
head(pc2_anti)
head(pc3_cor)
head(pc3_anti)

#  Combine PC1–PC3 gene sets

pc_genes <- unique(c(
  pc1_cor, pc1_anti,
  pc2_cor, pc2_anti,
  pc3_cor, pc3_anti
))
# Check number of unique genes extracted from PCA 
length(pc_genes)
# Inspect distribution of gene-PC correlations
summary(gene_pc_cor$PC1)
summary(gene_pc_cor$PC2)
summary(gene_pc_cor$PC3)
# Count number of genes passing correlation thresholds
sum(gene_pc_cor$PC1 > 0.2)
sum(gene_pc_cor$PC1 < -0.2)
# Confirm selected genes exist in filtered expression matrix
table(pc_genes %in% rownames(expr_var))
#  Identify Neural Progenitor Cells (NPCs)
# NPC markers from the paper:
#   PAX6, SOX2, VIM, HES1, PROM1
#
# Cells expressing ≥2 markers (>1 FPKM) are classified as NPCs

expr_cells <- t(expr_clean)

npc_markers <- c("PAX6","SOX2","VIM","HES1","PROM1")

npc_markers <- npc_markers[npc_markers %in% colnames(expr_cells)]

npc_score <- rowSums(expr_cells[, npc_markers, drop = FALSE] > 1)
# Inspect distribution of NPC marker scores 
table(npc_score)

# Identify cells classified as NPCs
npc_cells <- names(npc_score[npc_score >= 2])

length(npc_cells)
# Subset expression matrix to NPC cells only 
expr_npc <- expr_clean[, npc_cells, drop = FALSE]
# Inspect dimensions of NPC expression matrix 
dim(expr_npc)

# Filter genes for NPC PCA
#
#filtering rules:
#   expression >1 FPKM in >2 cells
#   variance >0.5

keep_expr_npc <- rowSums(expr_npc > 1) > 2

expr_expr_npc <- expr_npc[keep_expr_npc, , drop = FALSE]

gene_var_npc <- apply(expr_expr_npc, 1, var)

expr_var_npc <- expr_expr_npc[gene_var_npc > 0.5, , drop = FALSE]
# Inspect dimensions of filtered NPC gene matrix
dim(expr_var_npc)


# Prepare NPC PCA matrix
npc_input <- t(expr_var_npc)
# inspect dimenstions 
dim(npc_input)

# Remove problematic genes before PCA
#
# Remove genes that:
# contain NA or Inf values
# have zero variance
# This prevents PCA/SVD failures

# Identify genes suitable for PCA 
keep_good_npc <- apply(npc_input, 2, function(x){
  all(is.finite(x)) && sd(x) > 0
})
# Check number of genes retained for NPC PCA 
table(keep_good_npc)
# Build cleaned NPC PCA input matrix 
npc_input_clean <- npc_input[, keep_good_npc, drop = FALSE]
# Inspect cleaned NPC PCA matrix dimensions 
dim(npc_input_clean)


# Run PCA on NPC cells

npc_pca <- prcomp(
  npc_input_clean,
  center = TRUE,
  scale. = FALSE
)
# Inspect variance explained by NPC principal components 
summary(npc_pca)$importance[,1:5] 
# Extract NPC PCA genes
# Extract gene loadings from NPC PCA
npc_loadings <- npc_pca$rotation



# PC1 strongest genes
npc_pc1 <- rownames(npc_loadings)[
  order(abs(npc_loadings[,1]), decreasing = TRUE)
][1:50]



# PC2 strongest genes
npc_pc2 <- rownames(npc_loadings)[
  order(abs(npc_loadings[,2]), decreasing = TRUE)
][1:50]


# Combine PC1 and PC2 NPC gene lists 
npc_gene_list <- unique(c(npc_pc1, npc_pc2))
# Check number of unique NPC PCA genes 
length(npc_gene_list)


# Create combined gene set used for clustering analysis
#
# Combine genes from:
#   fetal PCA (PC1–PC3)
#   NPC PCA

combined_genes <- unique(c(
  pc1_cor, pc1_anti,
  pc2_cor, pc2_anti,
  pc3_cor, pc3_anti,
  npc_gene_list
))
# Inspect combined gene set size 
length(combined_genes)
head(npc_pc1)
head(npc_pc2)
# Build heatmap expression matrix
# Keep only genes present in the expression matrix

genes_present <- intersect(combined_genes, rownames(expr_clean))
# Check number of genes retained for heatmap 
length(genes_present)
# build heatmap matrix 
heatmap_mat <- t(expr_clean[genes_present, , drop = FALSE])
# inspect dimenstions 
dim(heatmap_mat)

# Remove problematic genes
#   Remove genes that would break scaling or clustering:
#   1. non-finite values
#   2. zero variance across cells
#
# Why:
#   scale() and correlation-based clustering fail or behave
#   poorly when genes are constant or contain invalid values.
# Identify genes suitable for clustering
keep_good_heatmap <- apply(heatmap_mat, 2, function(x) {
  all(is.finite(x)) && sd(x) > 0
})
# Check number of genes kept after filtering
table(keep_good_heatmap)
# Build cleaned heatmap matrix 
heatmap_mat <- heatmap_mat[, keep_good_heatmap, drop = FALSE]
# inspect 
dim(heatmap_mat)

# Remove marker genes so they appear only in the marker panel

marker_overlap <- intersect(gene_group_df2$gene, markers)
# Remove marker genes from gene grouping table
gene_group_df2 <- gene_group_df2[
  !gene_group_df2$gene %in% marker_overlap,
]
# Remove marker genes from heatmap expression matrix 
heatmap_log2 <- heatmap_log2[
  ,
  !colnames(heatmap_log2) %in% marker_overlap
]
# Compute log2-transformed expression matrix 
expr_log2 <- log2(expr_clean + 1)

# transpose to cells × genes
heatmap_log2 <- t(expr_log2)

# reorder cells according to clustering
heatmap_log2 <- heatmap_log2[rownames(heatmap_scaled2), ]

# reorder genes according to GO modules
main_mat <- heatmap_log2[, gene_order]

dim(main_mat)
# Compute gene order based on module assignments
gene_order <- order(gene_group_df2$group)
# Build final ordered heatmap matrix
main_mat <- heatmap_log2[, gene_order]
# Extract gene module labels for annotation
group_vec <- gene_group_df2$group[gene_order]
# Verify marker genes are excluded from main heatmap 
intersect(colnames(main_mat), markers)

#Hierarchical clustering of cells
#
#   distance = 1 - Pearson correlation
#   clustering = hierarchical clustering

cell_cor <- cor(t(heatmap_scaled), method = "pearson")
# Convert correlation matrix into clustering distance

cell_dist <- as.dist(1 - cell_cor)
# Perform hierarchical clustering of cells 
cell_hclust <- hclust(cell_dist, method = "average")

# inspect dendrogram
plot(cell_hclust, labels = FALSE, main = "Cell dendrogram")


#Inspect cut heights on dendrogram
#Visualize possible cut heights before assigning clusters.
# inspect the dendrogram and decide where the major
#   branches correspond to AP1, AP2, BP1, BP2, N1, N2, N3.
# Visualize clustering dendrogram 
plot(cell_hclust, labels = FALSE, main = "Cell dendrogram with suggested cut")
abline(h = 0.35, col = "red", lwd = 2)


# Assign cell clusters

cut_height <- 0.35

cell_clusters <- cutree(cell_hclust, k = 7)
# Inspect number of cells per cluster
table(cell_clusters)

# number of clusters obtained
length(unique(cell_clusters))

# Order cells by dendrogram
#   Use the dendrogram order for both the main heatmap and
#   the marker heatmap so all panels are aligned.

cell_order <- cell_hclust$order
# Reorder heatmap rows according to dendrogram order 
heatmap_ordered <- heatmap_scaled[cell_order, , drop = FALSE]

dim(heatmap_ordered)


# Prepare cell annotation

annotation_cells <- data.frame(
  Cluster = factor(cell_clusters[rownames(heatmap_ordered)])
)

rownames(annotation_cells) <- rownames(heatmap_ordered)

# Plot main heatmap

library(pheatmap)

pheatmap(
  t(heatmap_ordered),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  annotation_col = annotation_cells,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  border_color = NA
)



# Define marker genes for cell-type interpretation
markers <- c(
  "PAX6", "SOX2", "VIM", "HES1", "PROM1",
  "EOMES", "HES6", "INSM1", "ASPM",
  "MYT1L", "NEUROD6", "TBR1", "BCL11B", "DCX"
)
# Keep only markers present in the expression matrix 
markers_present <- intersect(markers, rownames(expr_clean))
# Inspect markers detected
markers_present
# Extract marker expression across cells
marker_mat <- t(expr_clean[markers_present, , drop = FALSE])

# Reorder marker matrix using the same cell order as heatmap
marker_mat <- marker_mat[rownames(heatmap_ordered), , drop = FALSE]

dim(marker_mat)
# Scale marker matrix for visualization
# Purpose: Standardize marker genes across cells so marker patterns
#are easily comparable across the ordered cells.

marker_scaled <- scale(marker_mat)

marker_scaled[marker_scaled > 3] <- 3
marker_scaled[marker_scaled < -3] <- -3

marker_scaled <- t(marker_scaled)


# Plot marker heatmap

pheatmap(
  marker_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  border_color = NA
)


# Save cell-cluster table

cluster_df <- data.frame(
  Cell = names(cell_clusters),
  Cluster = as.integer(cell_clusters)
)

head(cluster_df)
# Inspect number of cells per cluster
table(cluster_df$Cluster)

#Load bulk cortical zone reference data

path_ref_zones <- "/home/nhsas1/Critical Review/GSE38805_human_FPKM.txt"

zones_raw <- read.table(
  path_ref_zones,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

zones_13 <- zones_raw[, c("w13_VZ","w13_ISVZ","w13_OSVZ","w13_CP")]

colnames(zones_13) <- c("VZ","iSVZ","oSVZ","CP")

dim(zones_13)


#Convert Ensembl IDs to gene symbols
# Bulk dataset uses Ensembl IDs while our matrix 
# uses gene symbols. We convert Ensembl → SYMBOL.

library(org.Hs.eg.db)
library(AnnotationDbi)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(zones_13),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# remove genes without symbol
keep <- !is.na(gene_symbols)

zones_13 <- zones_13[keep, ]
gene_symbols <- gene_symbols[keep]

# remove duplicated symbols
unique_idx <- !duplicated(gene_symbols)

zones_13 <- zones_13[unique_idx, ]
gene_symbols <- gene_symbols[unique_idx]

rownames(zones_13) <- gene_symbols

dim(zones_13)

# Identify shared genes
# Correlation requires genes present in both datasets

shared_genes <- intersect(
  rownames(expr_clean),
  rownames(zones_13)
)

length(shared_genes)


# Restrict datasets to shared genes


expr_shared <- expr_clean[shared_genes, , drop = FALSE]
zones_shared <- zones_13[shared_genes, , drop = FALSE]

dim(expr_shared)
dim(zones_shared)



#Reorder cells according to clustering
# This ensures correlation results follow the same
# dendrogram order used in the heatmap.


expr_shared <- expr_shared[, rownames(heatmap_ordered), drop = FALSE]

dim(expr_shared)

#Remove problematic genes


# remove NA
keep_noNA <- complete.cases(zones_shared)

expr_shared  <- expr_shared[keep_noNA, ]
zones_shared <- zones_shared[keep_noNA, ]

# remove zero variance genes
var_expr <- apply(expr_shared, 1, var)
var_zone <- apply(zones_shared, 1, var)

keep_var <- var_expr > 0 & var_zone > 0

expr_shared  <- expr_shared[keep_var, ]
zones_shared <- zones_shared[keep_var, ]

dim(expr_shared)
dim(zones_shared)

#Compute cell vs zone correlations

zone_cor <- cor(expr_shared, zones_shared, method = "pearson")

dim(zone_cor)

head(zone_cor)


#Assign most correlated zone to each cell


zone_label <- colnames(zone_cor)[max.col(zone_cor)]

table(zone_label)


# Step 53 — Create zone annotation table

zone_annotation <- data.frame(
  Zone = factor(zone_label, levels = c("VZ","iSVZ","oSVZ","CP"))
)

rownames(zone_annotation) <- rownames(zone_cor)

head(zone_annotation)
table(cell_clusters, zone_label)

# Build gene group table from DAVID enrichment results

extract_david_genes <- function(file){
  
  lines <- readLines(file)
  
  gene_lines <- grep("GOTERM", lines, value = TRUE)
  
  split_lines <- strsplit(gene_lines, "\t")
  
  gene_strings <- sapply(split_lines, function(x) x[6])
  
  genes <- unique(unlist(strsplit(gene_strings, ", ")))
  
  return(genes)
}

# Load gene lists


pc1cor_genes  <- extract_david_genes("/home/nhsas1/RNA seq project /pc1cor.txt")
pc1anti_genes <- extract_david_genes("/home/nhsas1/RNA seq project /pc1anti.txt")

pc2cor_genes  <- extract_david_genes("/home/nhsas1/RNA seq project /PC2cor.txt")
pc2anti_genes <- extract_david_genes("/home/nhsas1/RNA seq project /PC2anti.txt")

pc2rganti_genes <- extract_david_genes("/home/nhsas1/RNA seq project /PC2rganti.txt")

pc3_genes <- extract_david_genes("/home/nhsas1/RNA seq project /PC3anti_and_PC2rgcor.txt")
pc4_genes <- extract_david_genes("/home/nhsas1/RNA seq project /PC4anti_and_PC1rgcor.txt")


# Define GO groups A–G 

gene_groups <- list(
  A = pc1cor_genes,
  B = pc1anti_genes,
  C = pc2cor_genes,
  D = pc2anti_genes,
  E = pc2rganti_genes,
  F = pc3_genes,
  G = pc4_genes
)


# Create gene-group annotation table


gene_group_df <- data.frame(
  gene = colnames(heatmap_scaled),
  group = "Z",
  stringsAsFactors = FALSE
)

for(g in names(gene_groups)){
  
  gene_group_df$group[
    gene_group_df$gene %in% gene_groups[[g]]
  ] <- g
  
}


# Heatmap


library(ComplexHeatmap)
library(circlize)
library(grid)


#Restrict genes to GO groups A–G


gene_group_df2 <- gene_group_df[
  gene_group_df$gene %in% colnames(heatmap_scaled),
]

gene_group_df2 <- gene_group_df2[
  match(colnames(heatmap_scaled), gene_group_df2$gene),
]

# keep only A–G modules
gene_group_df2$group[!(gene_group_df2$group %in% c("A","B","C","D","E","F","G"))] <- NA

keep_genes <- !is.na(gene_group_df2$group)

heatmap_scaled2 <- heatmap_scaled[, keep_genes]
gene_group_df2  <- gene_group_df2[keep_genes, ]


#Order genes by GO module


gene_group_df2$group <- factor(
  gene_group_df2$group,
  levels = c("A","B","C","D","E","F","G")
)

gene_order <- order(gene_group_df2$group)

main_mat <- heatmap_scaled2[, gene_order]

group_vec <- gene_group_df2$group[gene_order]
# Representative genes for each module (same as paper)
gene_labels <- c(
  A = "PAX6, GLI3,\nVIM, SOX2",
  B = "MYT1L, TBR1,\nNEUROD6, DCX",
  C = "NRCAM, FAT3,\nCDH6/7, SNCA",
  D = "NRP1, POU3F2,\nROBO2, ELMO1",
  E = "SEMA5A,\nDMD, PDPN",
  F = "MKI67,\nKIF11",
  G = "NEUROD4,\nEOMES"
)

bottom_annotation <- HeatmapAnnotation(
  genes = anno_block(
    labels = gene_labels,
    labels_gp = gpar(fontsize = 9)
  ),
  which = "column",
  height = unit(2.5, "cm")
)

# Prepare marker gene panel

marker_keep <- markers[markers %in% rownames(expr_clean)]

marker_mat <- t(
  expr_clean[marker_keep, rownames(main_mat), drop = FALSE]
)

# scale markers
marker_mat <- scale(marker_mat)

# clip extreme values
marker_mat[marker_mat > 3] <- 3
marker_mat[marker_mat < -3] <- -3

# Zone annotation (left bar)

names(zone_label) <- rownames(zone_cor)
zone_vec <- zone_label[rownames(main_mat)]
length(intersect(rownames(main_mat), names(zone_label)))
zone_annotation <- rowAnnotation(
  Zone = zone_vec,
  col = list(
    Zone = c(
      VZ   = "#FFD92F",
      iSVZ = "#66C2A5",
      oSVZ = "#FC8D62",
      CP   = "#8DA0CB"
    )
  )
)


#GO module labels (top annotation)


go_labels <- c(
  "Cell cycle\nForebrain development",
  "Neuron differentiation\n& projection",
  "Cell adhesion\nVesicle transport",
  "Neurogenesis\nCell migration",
  "Cell adhesion\nCell morphogenesis",
  "Cell cycle\nMitosis",
  "Neurogenesis\nCell morphogenesis"
)

top_annotation <- HeatmapAnnotation(
  
  GO = anno_block(
    labels = go_labels,
    labels_gp = gpar(fontsize = 8),
    labels_rot = 0,
    labels_just = "center",
    gp = gpar(fill = "grey90")
  ),
  
  height = unit(7, "cm")
)


# Main gene expression heatmap

# limit expression range
main_mat[main_mat > 6] <- 6
col_fun <- colorRamp2(
  c(0, 3, 6),
  c("navy", "white", "gold")
)
ht_main <- Heatmap(
  bottom_annotation = bottom_annotation , 
  main_mat,
  
  name = "Log2 FPKM",
  
  col = col_fun,
  
  cluster_rows = as.dendrogram(cell_hclust),
  cluster_columns = FALSE,
  
  column_split = factor(group_vec, levels = c("A","B","C","D","E","F","G")),
  
  column_gap = unit(18, "mm"),
  
  width = unit(ncol(main_mat) * 2.2, "mm"),
  
  show_row_names = FALSE,
  show_column_names = FALSE,
  
  left_annotation = zone_annotation,
  top_annotation = top_annotation
)


# Marker gene panel


ht_marker <- Heatmap(
  
  marker_mat,
  
  name = "Markers",
  
  col = colorRamp2(
    c(-3,0,3),
    c("#2166AC","white","#B2182B")
  ),
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  show_row_names = FALSE,
  show_column_names = TRUE,
  
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 10),
  
  width = unit(6.5, "cm")
  
)


# Draw final figure


draw(
  ht_main + ht_marker,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)


# 9. Save figure

pdf("Figure1C_reproduction.pdf", width = 20, height = 11)

draw(
  ht_main + ht_marker,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()
### final inspection 
ncol(expr_clean)
nrow(expr_clean)     # total genes
nrow(expr_expr)      # expressed genes
nrow(expr_var)       # variable genes
length(combined_genes)
length(shared_genes)

dim(expr_shared)
dim(zones_shared)
markers <- c(
  "PAX6","SOX2","VIM","HES1","PROM1",
  "EOMES","HES6","INSM1","ASPM",
  "MYT1L","NEUROD6","TBR1","BCL11B","DCX"
)

setdiff(markers, rownames(expr_clean))
table(cell_clusters)
summary(as.vector(log2(expr_clean + 1)))
summary(pca_res)$importance[,1:5]
library(cluster)

sil <- silhouette(cell_clusters, cell_dist)

mean(sil[,3])
summary(as.vector(zone_cor))
setdiff(rownames(expr_clean), rownames(zones_13))[1:20]
length(shared_genes) / nrow(expr_clean)
## 
comparison_summary <- data.frame(
  
  Metric = c(
    "Cells",
    "Total genes",
    "Expressed genes",
    "Variable genes",
    "PCA genes",
    "Shared genes with bulk"
  ),
  
  Your_Data = c(
    ncol(expr_clean),
    nrow(expr_clean),
    nrow(expr_expr),
    nrow(expr_var),
    length(combined_genes),
    length(shared_genes)
  )
  
)

print(comparison_summary)
paper_srr <- meta_table$SRR
my_cells <- colnames(expr_clean)
setdiff(paper_srr, my_cells)
setdiff(my_cells, paper_srr)
length(intersect(paper_srr, my_cells))
missing_srr <- setdiff(paper_srr, my_cells)

meta_table[meta_table$SRR %in% missing_srr, ]
endo_cell
missing_srr == endo_cell
paper_srr <- meta_table$SRR[meta_table$type == "fetal"]

my_cells <- colnames(expr_clean)

setdiff(paper_srr, my_cells)
length(paper_srr)
length(my_cells)
missing_srr <- setdiff(paper_srr, my_cells)

missing_srr
meta_table[meta_table$SRR %in% missing_srr, ]
cell_summary <- data.frame(
  
  Metric = c(
    "Fetal SRR runs in metadata",
    "Cells removed during filtering",
    "Cells used for analysis"
  ),
  
  Count = c(
    length(paper_srr),
    length(missing_srr),
    length(my_cells)
  )
  
)

cell_summary
removed_cells <- setdiff(paper_srr, my_cells)
removed_cells
