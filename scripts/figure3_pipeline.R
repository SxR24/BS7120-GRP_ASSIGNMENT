############################################################
# Figure 3 Gene Expression Analysis Pipeline
############################################################

# Load libraries
library(data.table)
library(pheatmap)

############################################################
# Set working directory
############################################################

setwd("~/figure3_work")

############################################################
# Load expression matrix
############################################################

expr <- fread("inputs/2026_UoL_prj_gene_matrix_f.tsv")

gene_names <- expr[[1]]
expr[[1]] <- NULL

expr_mat <- as.matrix(expr)
expr_mat <- apply(expr_mat, 2, as.numeric)

rownames(expr_mat) <- gene_names

############################################################
# Inspect matrix
############################################################

dim(expr_mat)
head(expr_mat[,1:5])

############################################################
# PCA analysis
############################################################

gene_var <- apply(expr_mat, 1, var, na.rm = TRUE)
gene_var[is.na(gene_var)] <- 0

expr_top <- expr_mat[
  gene_var > quantile(gene_var, 0.9, na.rm = TRUE),
  ,
  drop = FALSE
]

# remove genes with zero variance before PCA
expr_top <- expr_top[apply(expr_top, 1, var, na.rm = TRUE) > 0, , drop = FALSE]

pca <- prcomp(t(expr_top), scale.=TRUE)

png("figures/PCA_samples.png", width=1200, height=900)

plot(pca$x[,1], pca$x[,2],
     xlab="PC1",
     ylab="PC2",
     main="PCA of Samples",
     pch=19,
     col="darkblue")

dev.off()

############################################################
# Sample clustering
############################################################

dist_mat <- dist(t(expr_top))
hc <- hclust(dist_mat)

png("figures/sample_clustering.png", width=1200, height=900)

plot(hc, main="Sample Clustering")

dev.off()

############################################################
# Marker genes
############################################################

marker_genes <- c(
  "PAX6","SOX2","VIM","HES1","PROM1",
  "EOMES","HES6","INSM1","ASPM",
  "MYT1L","NEUROD6","TBR1","BCL11B","DCX"
)

writeLines(marker_genes,"outputs/marker_genes.tsv")

genes_present <- intersect(marker_genes, rownames(expr_mat))

writeLines(genes_present,"outputs/genes_heatmap.tsv")

############################################################
# Basic heatmap
############################################################

heatmap_mat <- expr_mat[genes_present,,drop=FALSE]

heatmap_scaled <- t(scale(t(heatmap_mat)))
heatmap_scaled[is.na(heatmap_scaled)] <- 0

png("figures/marker_heatmap.png", width=1400, height=1000)

pheatmap(
  heatmap_scaled,
  cluster_rows=TRUE,
  cluster_cols=TRUE,
  show_colnames=FALSE,
  fontsize_row=10,
  main="Marker Gene Heatmap"
)

dev.off()

############################################################
# Load metadata
############################################################

meta <- read.csv("inputs/SraRunTable (1).csv", stringsAsFactors=FALSE)

############################################################
# Match metadata to matrix
############################################################

common_ids <- intersect(meta$Run, colnames(expr_mat))

meta2 <- meta[match(common_ids, meta$Run),]

expr_mat2 <- expr_mat[,common_ids,drop=FALSE]

############################################################
# Annotation table
############################################################

annotation_col <- data.frame(
  stage = meta2$stage,
  tissue = meta2$tissue
)

rownames(annotation_col) <- meta2$Run

############################################################
# Annotated heatmap
############################################################

heatmap_mat2 <- expr_mat2[genes_present,,drop=FALSE]

heatmap_scaled2 <- t(scale(t(heatmap_mat2)))
heatmap_scaled2[is.na(heatmap_scaled2)] <- 0

png("figures/marker_heatmap_annotated.png", width=1800, height=1200)

pheatmap(
  heatmap_scaled2,
  annotation_col=annotation_col,
  cluster_rows=TRUE,
  cluster_cols=TRUE,
  show_colnames=FALSE,
  fontsize_row=10,
  main="Annotated Marker Gene Heatmap"
)

dev.off()

############################################################
# Save outputs
############################################################

write.csv(meta2,"outputs/metadata_aligned.csv",row.names=FALSE)

write.table(annotation_col,
            file="outputs/cell_metadata.tsv",
            sep="\t",
            quote=FALSE,
            col.names=NA)

############################################################
# Done
############################################################

print("Pipeline complete.")
