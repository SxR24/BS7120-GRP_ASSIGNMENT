# -----------------------------------------------------------------------------
# Reproduction of FigS3H (adapted for Seurat object)
# -----------------------------------------------------------------------------

library(Seurat)          
library(gplots)
library(RColorBrewer)


# Inputs
seurat_obj <- readRDS("seurat_preprocessed.rds")

expr <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
expr <- as.matrix(expr)

meta <- seurat_obj@meta.data


sample_clusters <- split(rownames(meta), meta$cluster)

# Marker gene sets

marker_sets <- list(
  `Cycling Cells` = c("ASPM", "UBE2C", "CDCA8", "CDC20", "CCNB1", "CCNB2",
                      "ETC2", "SGOL2", "TOP2A", "HMGB2"),
  Neuronal = c("SOX4", "DCX", "MLLT11", "KIF5A", "STMN2", "MYT1L",
               "RTN1", "NRXN1", "GAP43", "HMP19", "NEUROD6"),
  `Dorsal forebrain` = c("BCL11A", "NFIA", "NFIB", "C1orf61", "FABP7",
                         "CREB5", "IFI44L", "LHX2", "TFAP2C"),
  `Ventral forebrain` = c("IGDCC3", "OTX2", "COL2A1", "PRTG", "DLK1",
                          "CRABP1", "MGST1", "BCAT1", "MDK", "C11orf31"),
  `RSPO+` = c("RSPO2", "RSPO3", "WLS", "MAF", "NTRK2", "GJA1", "TPBG",
              "CTGF", "C1orf168", "RSPO1", "WNT8B"),
  Mesenchyme = c("LUM", "DCN", "COL3A1", "COL1A2", "COL1A1", "COL5A1",
                 "S100A10", "S100A11", "PITX2", "POSTN")
)

block_order <- names(marker_sets)


# Parse sample clusterss

sample_to_cluster <- stack(sample_clusters)
names(sample_to_cluster) <- c("Sample", "Cluster")
sample_to_cluster <- sample_to_cluster[sample_to_cluster$Sample %in% colnames(expr), ]

cluster_order <- names(sample_clusters)
ordered_samples <- c()
for (clust in cluster_order) {
  samples_in_clust <- intersect(sample_clusters[[clust]], colnames(expr))
  ordered_samples <- c(ordered_samples, samples_in_clust)
}
expr <- expr[, ordered_samples]

# Gene selection and matrix construction
selected_genes <- c()
block_labels <- c()

for (block in block_order) {
  present <- intersect(marker_sets[[block]], rownames(expr))
  if (length(present) > 0) {
    selected_genes <- c(selected_genes, present)
    block_labels <- c(block_labels, rep(block, length(present)))
  } else {
    message("Warning: No genes from block '", block, "' found in matrix.")
  }
}

mat <- as.matrix(expr[selected_genes, ])

# Row clustering within each block
row_order <- c()
for (block in block_order) {
  idx <- which(block_labels == block)
  block_mat <- mat[idx, , drop = FALSE]
  
  if (nrow(block_mat) > 1) {
    # Pearson correlation distance
    cor_mat <- tryCatch(
      cor(t(block_mat), method = "pearson", use = "pairwise.complete.obs"),
      error = function(e) NULL
    )
    if (!is.null(cor_mat)) {
      row_dist <- as.dist(1 - cor_mat)
      row_hclust <- hclust(row_dist, method = "average")
      ord <- row_hclust$order
      row_order <- c(row_order, idx[ord])
    } else {
      row_order <- c(row_order, idx)
    }
  } else {
    row_order <- c(row_order, idx)
  }
}

if (length(row_order) != length(selected_genes)) {
  warning("Row order length (", length(row_order), 
          ") does not match number of selected genes (", length(selected_genes), ").")
}

mat_final <- mat[row_order, ]
final_block_labels <- block_labels[row_order]

block_change <- which(final_block_labels[-1] != final_block_labels[-length(final_block_labels)])

# Column grouping and custom ordering 
group_definitions <- list(
  Group1 = c("C3", "C5", "C7"),
  Group2 = c("C1", "C2", "C4"),
  Group3 = c("C8"),
  Group4 = c("C6")
)

group_order <- c()
used_clusters <- c()
for (g in names(group_definitions)) {
  clusters <- group_definitions[[g]]
  new_clusters <- setdiff(clusters, used_clusters)
  if (length(new_clusters) < length(clusters)) {
    warning("Duplicate cluster(s) found in group definitions: ", 
            paste(setdiff(clusters, new_clusters), collapse = ", "),
            ". They will be placed only in their first group.")
  }
  group_order <- c(group_order, new_clusters)
  used_clusters <- c(used_clusters, new_clusters)
}

missing_clusters <- setdiff(names(sample_clusters), used_clusters)
if (length(missing_clusters) > 0) {
  warning("The following clusters were not included in any group and will be placed at the end: ",
          paste(missing_clusters, collapse = ", "))
  group_order <- c(group_order, missing_clusters)
}

new_col_order <- c()
for (clust in group_order) {
  samples <- sample_clusters[[clust]]
  samples <- intersect(samples, colnames(mat_final))
  if (length(samples) > 0) {
    new_col_order <- c(new_col_order, samples)
  }
}

all_cols <- colnames(mat_final)
if (!setequal(new_col_order, all_cols)) {
  missing_cols <- setdiff(all_cols, new_col_order)
  warning("Some samples were not placed: ", paste(missing_cols, collapse = ", "),
          ". They will be appended at the end.")
  new_col_order <- c(new_col_order, missing_cols)
}

mat_final <- mat_final[, new_col_order]

col_block_vec <- sapply(colnames(mat_final), function(s) {
  sample_to_cluster$Cluster[sample_to_cluster$Sample == s]
})

# Compute column separator positions
colsep_positions <- c()
cumulative_cols <- 0
for (g in names(group_definitions)) {
  clusters_in_group <- intersect(group_definitions[[g]], unique(col_block_vec))
  if (length(clusters_in_group) == 0) next
  
  group_sample_count <- sum(col_block_vec %in% clusters_in_group)
  cumulative_cols <- cumulative_cols + group_sample_count
  
  if (cumulative_cols < ncol(mat_final)) {
    colsep_positions <- c(colsep_positions, cumulative_cols)
  }
}

cat("\n--- Diagnostic output ---\n")
cat("Row separator positions:", block_change, "\n")
cat("Column separator positions (between groups):", colsep_positions, "\n")
cat("Number of genes per block:\n")
print(table(final_block_labels))
cat("Number of samples per cluster:\n")
print(table(col_block_vec))
cat("--------------------------\n\n")


# Colour definitions

cluster_colour_manual <- c(
  "C1"  = "#9A9CF8",
  "C2"  = "#5F88E3",
  "C3"  = "#C3352B",
  "C4"  = "#05E6EB",
  "C5"  = "#D2751C",
  "C6"  = "#F99BFB",
  "C7"  = "#A8197D",
  "C8"  = "#9DF9A8"
)

missing_cols <- setdiff(unique(col_block_vec), names(cluster_colour_manual))
if (length(missing_cols) > 0) {
  warning("The following clusters do not have a manual colour assigned and will be given random colours: ",
          paste(missing_cols, collapse = ", "))
  random_cols <- setNames(rainbow(length(missing_cols), s = 0.8, v = 0.9), missing_cols)
  cluster_colour_manual <- c(cluster_colour_manual, random_cols)
}
sample_cluster_colours <- cluster_colour_manual[unique(col_block_vec)]

block_colours <- setNames(
  brewer.pal(length(block_order), "Dark2"),
  block_order
)
row_side_cols <- block_colours[final_block_labels]
col_side_cols <- cluster_colour_manual[col_block_vec]

# Heatmap colour scale
expr_values <- as.vector(mat_final)
expr_limits <- quantile(expr_values, c(0.01, 0.99), na.rm = TRUE)
col_breaks <- seq(expr_limits[1], expr_limits[2], length.out = 101)
my_palette <- colorRampPalette(c("darkblue", "white", "red2"))(100)

row_names <- rownames(mat_final)
max_char <- max(nchar(row_names))
right_margin <- max(8, max_char * 0.18)

# 8. Draw heatmap

heatmap.2(mat_final,
          Rowv = FALSE,                     
          Colv = FALSE,                     
          dendrogram = "none",              
          scale = "none",
          col = my_palette,
          breaks = col_breaks,
          RowSideColors = row_side_cols,
          ColSideColors = col_side_cols,
          labRow = row_names,               
          labCol = rep("", ncol(mat_final)),
          srtRow = 0,                       
          adjRow = c(0, 0.5),               
          margins = c(8, right_margin),     
          key = TRUE,
          keysize = 1.2,
          density.info = "none",
          trace = "none",
          rowsep = block_change,            
          colsep = colsep_positions,        
          main = "Gene Expression Heatmap by organoid clusters"
)

# Add gene block labels
n_genes <- nrow(mat_final)
block_midpoints <- c()
unique_blocks <- unique(final_block_labels)
for (blk in unique_blocks) {
  idx <- which(final_block_labels == blk)
  mid <- floor(median(idx))
  block_midpoints <- c(block_midpoints, mid)
}

y_positions <- n_genes - block_midpoints + 1   

usr <- par("usr")
x_label <- usr[1] - (usr[2] - usr[1]) * 0.03

text(x = x_label, 
     y = y_positions, 
     labels = unique_blocks, 
     col = block_colours[unique_blocks], 
     font = 2,          
     cex = 0.9, 
     adj = c(1, 0.5),   
     xpd = NA)          

# Add legends
natural_sort <- function(x) {
  x[order(as.numeric(gsub("[^0-9]", "", x)))]
}

sorted_clusters <- natural_sort(unique(col_block_vec))
sample_cluster_colours_sorted <- cluster_colour_manual[sorted_clusters]

par(lend = 1)
legend("bottomleft",
       legend = names(block_colours),
       fill = block_colours,
       border = NA,
       bty = "n",
       cex = 0.8,
       title = "Gene Blocks",
       inset = c(-0.05, -0.1),
       xpd = TRUE)

legend("left",
       legend = names(sample_cluster_colours_sorted),
       fill = sample_cluster_colours_sorted,
       border = NA,
       bty = "n",
       cex = 0.8,
       title = "Sample Clusters",
       inset = c(-0.05, -0.1),
       xpd = TRUE)