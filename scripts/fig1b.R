
# FIGURE 1B REPLICATION , Camp et al.,
# load the packages
suppressPackageStartupMessages({
  library(pheatmap)
  library(gridExtra)
})
# Stores file locations in variables to separate configuration
path_sc_master <- "/home/nhsas1/Critical Review/GSE75140_hOrg.fetal.master.data.frame.txt"
path_ref_zones <- "/home/nhsas1/Critical Review/GSE38805_human_FPKM.txt"
path_ref_cells <- "/home/nhsas1/Critical Review/GSE65000_hsa_fpkm_matrix.txt"

# Load the log2 FPKM single cell matrix from GEO (GSE75140) which contains data ofo both organoids and fetal cells ( Cells x Genes )

sc_raw <- read.table(path_sc_master,
                     header = TRUE,
                     row.names = 1,
                     sep = "\t",
                     check.names = FALSE)

# Keep only fetal cells 
sc_fetal <- sc_raw[grep("fetal", rownames(sc_raw)), ]

# Keep only 12–13 wpc 
sc_fetal <- sc_fetal[grep("12wpc|13wpc", rownames(sc_fetal)), ]

# Check the retained cells to verify correct developmental-stage filtering 
cat("Cells after fetal 12–13 filtering:", nrow(sc_fetal), "\n")

# Convert to numeric matrix and transpose → Genes as rows x Cells as cols 
sc_matrix <- t(data.matrix(sc_fetal))
mode(sc_matrix) <- "numeric"

# Remove endothelial (PECAM1) and five interneuron cells (GAD1,ERBB4,DLX1,DLX2,DLX5,DLX6)
# Store cell IDs 
cells <- colnames(sc_matrix)

#Helper function to safely extract expression values for a given gene , checks whether the gene symbol exists in the sc matrix
# before attempting row indexing, preventing subscript errors if a marker
# gene is absent. Returns the gene's expression vector across all cells or NULL if the gene is not present

get_gene <- function(gene) {
  if (gene %in% rownames(sc_matrix)) {
    return(sc_matrix[gene, ])
  } else {
    return(NULL)
  }
}

# Extract PECAM1 expression across all cells
pecam1 <- get_gene("PECAM1")
#Create a variable to store which cell will be removed
endothelial_cell <- character(0)
#Check whether PECAM1 exists and is expressed in at least one cell , 
#then pick “most PECAM1-high” cell and store its ID
if (!is.null(pecam1) && any(pecam1 > 0)){
  endothelial_cell <- names(which.max(pecam1))
}

# ---- Interneuron markers
# Make a list of interneurons marker genes
markers <- c("GAD1","ERBB4","DLX1","DLX2","DLX5","DLX6")
# Create a score vector: one score per cell, starting at 0
score <- setNames(rep(0, length(cells)), cells)
# Iterate over each interneuron marker gene. For each marker, retrieve its expression across all cells 
#and add it to the cumulative score if the gene is present.
for (m in markers) {
  gene_vec <- get_gene(m)
  if (!is.null(gene_vec)) {
    score <- score + gene_vec
  }
}
# Initialize container for identified interneuron cell IDs.
# Starts empty in case no cells show positive marker expression
interneurons <- character(0)
# Proceed only if at least one cell has a positive cumulative interneuron marker score
# then rank cells in descending order and select top 5 
if (any(score > 0)) {
  ranked <- names(sort(score, decreasing = TRUE))
  interneurons <- head(ranked[score[ranked] > 0], 5)
}
# Remove the identified cells
cells_to_remove <- unique(c(endothelial_cell, interneurons))
# Verify the total removed is 6 
cat("Cells removed:", length(cells_to_remove), "\n")
sc_matrix <- sc_matrix[, !colnames(sc_matrix) %in% cells_to_remove]
# Verify the remaining cells are 220 
cat("Cells remaining:", ncol(sc_matrix), "\n")


#Load cortical zone bulk reference
zones_raw <- read.table(path_ref_zones,
                        header = TRUE,
                        row.names = 1,
                        sep = "\t",
                        check.names = FALSE)
# Extract the wpc 13 only
zones_13 <- zones_raw[, c("w13_VZ","w13_ISVZ","w13_OSVZ","w13_CP")]
# Renames columns
colnames(zones_13) <- c("VZ","iSVZ","oSVZ","CP")


# Load purified cell-type reference

cells_raw <- read.table(path_ref_cells,
                        header = TRUE,
                        row.names = 1,
                        sep = "\t",
                        check.names = FALSE)
# Select cells cols 
aRG_cols <- grep("^Hsa_aRG", colnames(cells_raw))
bRG_cols <- grep("^Hsa_bRG", colnames(cells_raw))
N_cols   <- grep("^Hsa_N_",  colnames(cells_raw))
# Compute the mean of cells replicates into one profile per type
cells_mean <- data.frame(
  aRG = rowMeans(cells_raw[, aRG_cols]),
  bRG = rowMeans(cells_raw[, bRG_cols]),
  N   = rowMeans(cells_raw[, N_cols])
)


# Map Ensembl IDs → gene symbols
gene_map <- data.frame(
  ensembl = rownames(cells_raw),
  symbol  = cells_raw$external_gene_id,
  stringsAsFactors = FALSE
)
# Removes missing symbols and duplicated Ensembl IDs
gene_map <- gene_map[!is.na(gene_map$symbol) & gene_map$symbol != "", ]
gene_map <- gene_map[!duplicated(gene_map$ensembl), ]

map_to_symbol <- function(mat) {
  
  keep <- rownames(mat) %in% gene_map$ensembl
  mat  <- mat[keep, ]
  
  symbols <- gene_map$symbol[match(rownames(mat), gene_map$ensembl)]
  
  df <- as.data.frame(mat)
  df$symbol <- symbols
  
  collapsed <- aggregate(. ~ symbol, data = df, FUN = mean)
  
  rownames(collapsed) <- collapsed$symbol
  collapsed$symbol <- NULL
  
  as.matrix(collapsed)
}

zones_symbol <- map_to_symbol(zones_13)
cells_symbol <- map_to_symbol(cells_mean)



# Intersect genes across all datasets

common_genes <- Reduce(intersect,
                       list(rownames(sc_matrix),
                            rownames(zones_symbol),
                            rownames(cells_symbol))) 
# Verify the shared genes (9699)

cat("Shared genes:", length(common_genes), "\n")

sc_matrix     <- sc_matrix[common_genes, ]
zones_symbol  <- zones_symbol[common_genes, ]
cells_symbol  <- cells_symbol[common_genes, ]


# Spearman correlation (cell vs references)

num_cells <- ncol(sc_matrix)

zone_cor  <- matrix(NA, nrow = num_cells, ncol = 4)
type_cor  <- matrix(NA, nrow = num_cells, ncol = 3)

colnames(zone_cor) <- colnames(zones_symbol)
colnames(type_cor) <- colnames(cells_symbol)
rownames(zone_cor) <- colnames(sc_matrix)
rownames(type_cor) <- colnames(sc_matrix)

for (i in 1:num_cells) {
  
  cell_vector <- sc_matrix[, i]
  
  for (j in 1:4) {
    zone_cor[i,j] <- cor(cell_vector,
                         zones_symbol[, j],
                         method = "spearman")
  }
  
  for (j in 1:3) {
    type_cor[i,j] <- cor(cell_vector,
                         cells_symbol[, j],
                         method = "spearman")
  }
}


# Per-cell Z-score normalization 
zone_z <- t(scale(t(zone_cor)))
type_z <- t(scale(t(type_cor)))



# Cluster cells using zone panel (Pearson correlation distance) 
distance_matrix <- as.dist(1 - cor(t(zone_z)))
hc <- hclust(distance_matrix, method = "complete")

row_order <- hc$order
# Apply identical row order to both panels
zone_z  <- zone_z[row_order, ]
type_z  <- type_z[row_order, ]


# Row annotation: wpc12 vs wpc13

cell_ids <- rownames(zone_z)

annotation_row <- data.frame(
  wpc = ifelse(grepl("12wpc", cell_ids), "wpc12", "wpc13"),
  row.names = cell_ids
)



# Build the heatmap

z_lim <- 2
breaks <- seq(-z_lim, z_lim, length.out = 101)
colors <- colorRampPalette(c("#2166ac","white","#b2182b"))(100)


hm_zone <- pheatmap(zone_z,
                    cluster_rows = FALSE,
                    cluster_cols = FALSE,
                    show_rownames = FALSE,
                    annotation_row = annotation_row,
                    color = colors,
                    breaks = breaks,
                    border_color = NA,
                    silent = TRUE)

hm_type <- pheatmap(type_z,
                    cluster_rows = FALSE,
                    cluster_cols = FALSE,
                    show_rownames = FALSE,
                    annotation_row = annotation_row,
                    color = colors,
                    breaks = breaks,
                    border_color = NA,
                    silent = TRUE)
# Arrange two panels 
grid.arrange(hm_zone$gtable,
             hm_type$gtable,
             ncol = 2)

# save image 
save.image("Fig1B_workspace.RData")
