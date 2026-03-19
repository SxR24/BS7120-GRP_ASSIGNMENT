# app.R
# ------------------------------------------------------------
# GSE75140 Shiny Explorer (FALLBACK VERSION)
# - No plotly
# - No Seurat
# - Uses only precomputed embedding RDS files
# - Works with shiny + ggplot2 + dplyr + DT + bslib
# ------------------------------------------------------------

library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(dplyr)

# -----------------------------
# Files
# -----------------------------
TSNE_RDS <- "emb_tsne.rds"
UMAP_RDS <- "emb_umap.rds"
PCA_RDS  <- "emb_pca.rds"

# -----------------------------
# Helpers
# -----------------------------
make_palette <- function(levels_vec) {
  cols <- grDevices::hcl.colors(length(levels_vec), palette = "Dark 3")
  names(cols) <- levels_vec
  cols
}

apply_alignment <- function(df, xcol, ycol, swap_xy = FALSE, flip_x = FALSE, flip_y = FALSE) {
  x <- df[[xcol]]
  y <- df[[ycol]]
  
  if (isTRUE(swap_xy)) {
    tmp <- x
    x <- y
    y <- tmp
  }
  
  if (isTRUE(flip_x)) x <- -x
  if (isTRUE(flip_y)) y <- -y
  
  df$X_plot <- x
  df$Y_plot <- y
  df
}

cluster_centroids <- function(df, cluster_col = "cluster") {
  df %>%
    group_by(.data[[cluster_col]]) %>%
    summarise(
      Xc = median(X_plot, na.rm = TRUE),
      Yc = median(Y_plot, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(cluster = 1)
}

pick_default_axes <- function(df) {
  num_cols <- names(df)[sapply(df, is.numeric)]
  preferred <- num_cols[grepl("^(tSNE|UMAP|PC)_", num_cols)]
  
  if (length(preferred) >= 2) {
    list(
      choices = preferred,
      x = preferred[1],
      y = preferred[2]
    )
  } else {
    list(
      choices = num_cols,
      x = num_cols[1],
      y = num_cols[2]
    )
  }
}

# -----------------------------
# Load precomputed embeddings
# -----------------------------
emb_tsne <- readRDS(TSNE_RDS)
emb_umap <- readRDS(UMAP_RDS)
emb_pca  <- readRDS(PCA_RDS)

# Basic checks
required_cols <- c("cell_id", "cluster")
for (nm in c("emb_tsne", "emb_umap", "emb_pca")) {
  obj <- get(nm)
  miss <- setdiff(required_cols, colnames(obj))
  if (length(miss) > 0) {
    stop(paste0(nm, " is missing required columns: ", paste(miss, collapse = ", ")))
  }
}

# Make sure species exists even if absent
if (!("species" %in% colnames(emb_tsne))) emb_tsne$species <- NA_character_
if (!("species" %in% colnames(emb_umap))) emb_umap$species <- NA_character_
if (!("species" %in% colnames(emb_pca)))  emb_pca$species  <- NA_character_

cluster_levels <- sort(unique(as.character(emb_tsne$cluster)))
CLUST_COLORS <- make_palette(cluster_levels)

N_CELLS <- nrow(emb_tsne)
N_CLUST <- length(cluster_levels)

# -----------------------------
# Theme
# -----------------------------
theme_pretty <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  base_font = "Segoe UI",
  heading_font = "Segoe UI Semibold",
  code_font = "Consolas",
  primary = "#4f46e5",
  secondary = "#0ea5e9",
  success = "#22c55e",
  danger = "#ef4444",
  warning = "#f59e0b",
  info = "#06b6d4",
  base_border_radius = "14px"
)

# -----------------------------
# UI
# -----------------------------
ui <- navbarPage(
  title = tagList(
    tags$span(style = "font-weight:700; letter-spacing:0.2px;", "GSE75140 Explorer"),
    tags$span(style = "opacity:0.7; font-weight:600; margin-left:8px; font-size:0.92rem;", "• fallback mode")
  ),
  theme = theme_pretty,
  
  header = tags$head(
    tags$style(HTML("
      :root{
        --glass-bg: rgba(255,255,255,0.68);
        --glass-brd: rgba(15,23,42,0.10);
        --shadow: 0 10px 30px rgba(2,6,23,0.10);
        --shadow2: 0 6px 18px rgba(2,6,23,0.10);
        --nav-bg: rgba(255,255,255,0.92);
        --nav-text: #0f172a;
        --nav-muted: rgba(15,23,42,0.62);
        --nav-active: #16a34a;
        --nav-border: rgba(15,23,42,0.10);
      }

      html, body { height: 100%; }
      body {
        overflow: hidden;
        background:
          radial-gradient(1100px 520px at 10% 0%, rgba(79,70,229,0.18), transparent 55%),
          radial-gradient(900px 520px at 95% 10%, rgba(14,165,233,0.18), transparent 55%),
          linear-gradient(180deg, #f8fafc, #f1f5f9);
        color: #0f172a;
      }

      .navbar, .navbar .container-fluid {
        background: var(--nav-bg) !important;
        border-bottom: 1px solid var(--nav-border) !important;
        backdrop-filter: blur(10px);
      }

      .navbar .navbar-brand,
      .navbar .navbar-brand:hover,
      .navbar .navbar-brand:focus {
        color: var(--nav-text) !important;
        font-weight: 800;
      }

      .navbar-nav .nav-link,
      .navbar-nav .nav-link:visited {
        color: var(--nav-muted) !important;
        font-weight: 700;
      }

      .navbar-nav .nav-link:hover,
      .navbar-nav .nav-link:focus {
        color: var(--nav-text) !important;
      }

      .navbar-nav .nav-link.active,
      .navbar-nav .show > .nav-link {
        color: var(--nav-active) !important;
        border-bottom: 2px solid var(--nav-active);
        padding-bottom: 10px;
      }

      .sidebar-scroll {
        position: sticky;
        top: 72px;
        height: calc(100vh - 90px);
        overflow-y: auto;
        padding-right: 10px;
      }

      .bslib-sidebar {
        background: var(--glass-bg) !important;
        border: 1px solid var(--glass-brd) !important;
        box-shadow: var(--shadow2);
      }

      .main-fixed {
        height: calc(100vh - 90px);
        overflow: hidden;
        display: flex;
        align-items: center;
        justify-content: center;
        padding: 10px;
      }

      .card {
        border: 1px solid rgba(15,23,42,0.10) !important;
        box-shadow: var(--shadow);
        background: rgba(255,255,255,0.78) !important;
        backdrop-filter: blur(10px);
      }

      .card-header {
        background: transparent !important;
        border-bottom: 1px solid rgba(15,23,42,0.08) !important;
        font-weight: 800;
        color: #0f172a;
      }

      .card-fullheight { height: 100%; width: 100%; }

      .chiprow { display:flex; gap:10px; flex-wrap:wrap; margin-bottom:10px; }
      .chip {
        padding: 8px 10px;
        border-radius: 999px;
        border: 1px solid rgba(15,23,42,0.10);
        background: rgba(255,255,255,0.65);
        font-weight: 700;
        font-size: 0.9rem;
        box-shadow: 0 4px 10px rgba(2,6,23,0.06);
        color: #0f172a;
      }
      .chip b { font-weight: 900; }
    "))
  ),
  
  tabPanel(
    "Clusters",
    layout_sidebar(
      sidebar = sidebar(
        width = 430,
        div(
          class = "sidebar-scroll",
          
          div(
            class = "chiprow",
            div(class = "chip", HTML(paste0("Cells: <b>", N_CELLS, "</b>"))),
            div(class = "chip", HTML(paste0("Clusters: <b>", N_CLUST, "</b>")))
          ),
          
          h4("Reduction"),
          selectInput(
            "reduction",
            "Reduction",
            choices = c("tSNE" = "tsne", "UMAP" = "umap", "PCA" = "pca"),
            selected = "tsne"
          ),
          uiOutput("axis_selectors"),
          
          hr(),
          h4("Alignment"),
          checkboxInput("swap_xy", "Swap X/Y", value = FALSE),
          checkboxInput("flip_x", "Flip X", value = FALSE),
          checkboxInput("flip_y", "Flip Y", value = FALSE),
          
          hr(),
          h4("Style"),
          sliderInput("pt_size", "Point size", min = 0.5, max = 4, value = 1.2, step = 0.1),
          sliderInput("pt_alpha", "Point alpha", min = 0.1, max = 1.0, value = 0.85, step = 0.05),
          checkboxInput("show_labels", "Show cluster labels", value = TRUE),
          
          hr(),
          h4("Run summary"),
          verbatimTextOutput("run_summary")
        )
      ),
      
      div(
        class = "main-fixed",
        card(
          class = "card-fullheight",
          card_header("Cluster map"),
          plotOutput("cluster_plot", height = "calc(100vh - 170px)")
        )
      )
    )
  ),
  
  tabPanel(
    "Summary",
    layout_sidebar(
      sidebar = sidebar(
        width = 430,
        div(
          class = "sidebar-scroll",
          h4("Summary table"),
          selectInput(
            "summary_group",
            "Group by",
            choices = c("Cluster" = "cluster", "Species" = "species"),
            selected = "cluster"
          )
        )
      ),
      
      div(
        class = "main-fixed",
        card(
          class = "card-fullheight",
          card_header("Cell summary"),
          DTOutput("summary_table")
        )
      )
    )
  )
)

# -----------------------------
# Server
# -----------------------------
server <- function(input, output, session) {
  
  output$run_summary <- renderText({
    paste0(
      "Mode: fallback / no plotly / no Seurat\n",
      "Cells loaded: ", N_CELLS, "\n",
      "Clusters loaded: ", N_CLUST, "\n",
      "Embeddings available: tSNE, UMAP, PCA\n",
      "Gene view: disabled in fallback mode"
    )
  })
  
  embedding_df <- reactive({
    switch(
      input$reduction,
      "tsne" = emb_tsne,
      "umap" = emb_umap,
      "pca"  = emb_pca
    )
  })
  
  output$axis_selectors <- renderUI({
    df <- embedding_df()
    ax <- pick_default_axes(df)
    
    tagList(
      selectInput("x_axis", "X axis", choices = ax$choices, selected = ax$x),
      selectInput("y_axis", "Y axis", choices = ax$choices, selected = ax$y)
    )
  })
  
  plot_df <- reactive({
    req(input$x_axis, input$y_axis)
    
    df <- embedding_df()
    df$cluster <- factor(as.character(df$cluster), levels = names(CLUST_COLORS))
    
    apply_alignment(
      df,
      input$x_axis,
      input$y_axis,
      input$swap_xy,
      input$flip_x,
      input$flip_y
    )
  })
  
  output$cluster_plot <- renderPlot({
    df <- plot_df()
    
    p <- ggplot(df, aes(x = X_plot, y = Y_plot, color = cluster)) +
      geom_point(size = input$pt_size, alpha = input$pt_alpha) +
      scale_color_manual(values = CLUST_COLORS) +
      labs(
        title = paste0("Cluster map: ", toupper(input$reduction)),
        x = input$x_axis,
        y = input$y_axis,
        color = "Cluster"
      ) +
      theme_classic(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "right"
      )
    
    if (isTRUE(input$show_labels)) {
      cents <- cluster_centroids(df, "cluster")
      cents$label <- as.character(cents$cluster)
      
      p <- p +
        geom_text(
          data = cents,
          aes(x = Xc, y = Yc, label = label),
          inherit.aes = FALSE,
          color = "black",
          size = 5,
          fontface = "bold"
        )
    }
    
    p
  })
  
  output$summary_table <- renderDT({
    df <- emb_tsne
    
    if (!("species" %in% colnames(df))) {
      df$species <- NA_character_
    }
    
    if (input$summary_group == "cluster") {
      out <- df %>%
        mutate(cluster = as.character(cluster)) %>%
        count(cluster, name = "n_cells") %>%
        arrange(cluster)
    } else {
      out <- df %>%
        mutate(species = ifelse(is.na(species) | species == "", "NA", species)) %>%
        count(species, name = "n_cells") %>%
        arrange(species)
    }
    
    datatable(out, options = list(pageLength = 15, scrollX = TRUE))
  })
}

shinyApp(ui, server)