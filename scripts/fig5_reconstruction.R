setwd("~/figure5_work")

library(igraph)

# -----------------------------
# 1. Input data
# -----------------------------
seed_genes <- c(
  "PROM1","SOX2","PAX6","GLI3","INSM1","HES6",
  "SOX11","SOX4","NEUROD6","MYT1L","TBR1","MEF2C","EOMES"
)

edge_df <- read.csv("results/fig5_gene_correlation_edges_big.csv", stringsAsFactors = FALSE)
screen_df <- read.csv("data/sd04_correlation_screens_clean.csv", stringsAsFactors = FALSE)

# Clean category names
screen_df$set_name <- trimws(screen_df$set_name)
screen_df$set_name[screen_df$set_name == "hcondel"] <- "hCondel"
screen_df$set_name[screen_df$set_name == "neanderthal"] <- "modHum"

# -----------------------------
# 2. Build main network
# -----------------------------
g_big <- graph_from_data_frame(edge_df, directed = FALSE)
E(g_big)$weight <- edge_df$cor

comp <- components(g_big)
largest_comp_id <- which.max(comp$csize)
g_main <- induced_subgraph(g_big, vids = V(g_big)[comp$membership == largest_comp_id])

V(g_main)$is_seed <- V(g_main)$name %in% seed_genes

set.seed(123)
lay_main <- layout_with_fr(g_main)

# -----------------------------
# 3. Left panel
# -----------------------------
png("results/fig5_left_panel_BEST.png", width = 1800, height = 1600)

plot(
  g_main,
  layout = lay_main,
  vertex.size = ifelse(V(g_main)$is_seed, 5, 1.5),
  vertex.label = NA,
  vertex.color = ifelse(V(g_main)$is_seed, "red", "black"),
  vertex.frame.color = NA,
  edge.width = 0.05,
  edge.color = rgb(0.7, 0.7, 0.7, 0.10),
  margin = 0.03
)

seed_ids <- which(V(g_main)$is_seed)
text(
  lay_main[seed_ids, 1],
  lay_main[seed_ids, 2],
  labels = V(g_main)$name[seed_ids],
  col = "red",
  cex = 0.75,
  pos = 3
)

dev.off()

# -----------------------------
# 4. Right panel overlays
# -----------------------------
screen_net <- screen_df[screen_df$gene %in% V(g_main)$name, ]

make_overlay_colors <- function(set_name, highlight_col) {
  genes_use <- screen_net$gene[screen_net$set_name == set_name]
  cols <- rep("black", vcount(g_main))
  sizes <- rep(1.2, vcount(g_main))
  idx <- which(V(g_main)$name %in% genes_use)
  cols[idx] <- highlight_col
  sizes[idx] <- 3.6
  list(cols = cols, sizes = sizes)
}

modHum_info  <- make_overlay_colors("modHum",  "green3")
OMIM_info    <- make_overlay_colors("OMIM",    "turquoise3")
hCondel_info <- make_overlay_colors("hCondel", "royalblue3")
haDHS_info   <- make_overlay_colors("haDHS",   "darkorange2")

png("results/fig5_right_panels_BEST.png", width = 1800, height = 1800)
par(mfrow = c(2,2), mar = c(0.5,0.5,2,0.5))

plot(
  g_main, layout = lay_main,
  vertex.size = modHum_info$sizes,
  vertex.label = NA,
  vertex.color = modHum_info$cols,
  vertex.frame.color = NA,
  edge.width = 0.08,
  edge.color = rgb(0.7,0.7,0.7,0.12),
  main = "i  modHum"
)

plot(
  g_main, layout = lay_main,
  vertex.size = OMIM_info$sizes,
  vertex.label = NA,
  vertex.color = OMIM_info$cols,
  vertex.frame.color = NA,
  edge.width = 0.08,
  edge.color = rgb(0.7,0.7,0.7,0.12),
  main = "ii  OMIM"
)

plot(
  g_main, layout = lay_main,
  vertex.size = hCondel_info$sizes,
  vertex.label = NA,
  vertex.color = hCondel_info$cols,
  vertex.frame.color = NA,
  edge.width = 0.08,
  edge.color = rgb(0.7,0.7,0.7,0.12),
  main = "iii  hCondel"
)

plot(
  g_main, layout = lay_main,
  vertex.size = haDHS_info$sizes,
  vertex.label = NA,
  vertex.color = haDHS_info$cols,
  vertex.frame.color = NA,
  edge.width = 0.08,
  edge.color = rgb(0.7,0.7,0.7,0.12),
  main = "iv  haDHS"
)

dev.off()
