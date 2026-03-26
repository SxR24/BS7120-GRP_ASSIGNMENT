library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(dplyr)


# Precomputed files

TSNE_RDS <- "emb_tsne.rds"
UMAP_RDS <- "emb_umap.rds"
PCA_RDS  <- "emb_pca.rds"

# Helpers

make_palette <- function(levels_vec) {
  base_cols <- c(
    "#0B8F66", "#1AA6A8", "#0F766E", "#14B8A6", "#0EA5E9",
    "#2563EB", "#4F46E5", "#7C3AED", "#DC2626", "#F59E0B",
    "#84CC16", "#475569", "#E11D48", "#0891B2", "#9333EA"
  )
  cols <- rep(base_cols, length.out = length(levels_vec))
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
    if (length(num_cols) < 2) {
      stop("Embedding file must contain at least two numeric columns for plotting.")
    }
    list(
      choices = num_cols,
      x = num_cols[1],
      y = num_cols[2]
    )
  }
}

safe_species <- function(df) {
  if (!("species" %in% colnames(df))) {
    df$species <- "NA"
  }
  df$species[is.na(df$species) | df$species == ""] <- "NA"
  df
}


# Load data

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

emb_tsne <- safe_species(emb_tsne)
emb_umap <- safe_species(emb_umap)
emb_pca  <- safe_species(emb_pca)

cluster_levels <- sort(unique(as.character(c(
  emb_tsne$cluster,
  emb_umap$cluster,
  emb_pca$cluster
))))
CLUST_COLORS <- make_palette(cluster_levels)

N_CELLS <- nrow(emb_tsne)
N_CLUST <- length(cluster_levels)
species_exists <- "species" %in% colnames(emb_tsne) ||
                  "species" %in% colnames(emb_umap) ||
                  "species" %in% colnames(emb_pca)


# Theme

theme_pretty <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  base_font = font_google("Inter"),
  heading_font = font_google("Inter"),
  code_font = "Consolas",
  primary = "#0B8F66",
  secondary = "#1AA6A8",
  success = "#16a34a",
  info = "#0891b2",
  warning = "#f59e0b",
  danger = "#dc2626",
  bg = "#F7FBF9",
  fg = "#0F172A",
  base_border_radius = "16px",
  btn_border_radius = "12px",
  input_border_radius = "12px"
)


# UI

ui <- navbarPage(
  id = "main_tabs",
  title = div(
    style = "display:flex; align-items:center; gap:12px;",
    tags$div(
      style = "
        width:14px; height:14px; border-radius:50%;
        background: linear-gradient(135deg, #DFF7ED, #7ADDBB);
        box-shadow: 0 0 14px rgba(255,255,255,0.35);
      "
    ),
    div(
      tags$span(
        style = "font-weight:800; letter-spacing:0.25px; color:white;",
        "Group A Dashboard"
      ),
      tags$span(
        style = "margin-left:12px; color:rgba(255,255,255,0.85); font-weight:500; font-size:0.92rem;",
        ""
      )
    )
  ),
  theme = theme_pretty,
  
  header = tags$head(
    tags$style(HTML("
      :root{
        --ebi-green: #0B8F66;
        --ebi-green-dark: #087556;
        --ebi-green-soft: #EAF7F2;
        --ebi-teal: #1AA6A8;
        --ink: #0F172A;
        --muted: #64748B;
        --line: #DCE7E5;
        --panel: rgba(255,255,255,0.92);
        --shadow-soft: 0 10px 30px rgba(15, 23, 42, 0.07);
        --shadow-card: 0 8px 24px rgba(15, 23, 42, 0.08);
        --ring: rgba(26,166,168,0.16);
      }

      html, body { height: 100%; }
      body {
        overflow: hidden;
        color: var(--ink);
        background:
          radial-gradient(900px 500px at 0% 0%, rgba(11,143,102,0.08), transparent 60%),
          radial-gradient(800px 420px at 100% 0%, rgba(26,166,168,0.08), transparent 60%),
          linear-gradient(180deg, #F8FCFB 0%, #F4F9F7 100%);
        font-family: 'Inter', sans-serif;
      }

      .navbar, .navbar .container-fluid {
        background: linear-gradient(90deg, var(--ebi-green-dark), var(--ebi-green)) !important;
        border-bottom: none !important;
        min-height: 40px;
        box-shadow: 0 4px 18px rgba(11,143,102,0.20);
      }

      .navbar-brand,
      .navbar-brand:hover,
      .navbar-brand:focus {
        color: white !important;
        font-weight: 800;
      }

      .navbar-nav {
        gap: 18px !important;
        margin-left: 400px !important;
      }

      .navbar-nav .nav-item {
        margin-right: 0px !important;
      }

      .navbar-nav .nav-link,
      .navbar-nav .nav-link:visited {
        color: rgba(255,255,255,0.88) !important;
        font-weight: 700;
        letter-spacing: 0.2px;
        padding-left: 14px !important;
        padding-right: 14px !important;
        padding-top: 12px !important;
        padding-bottom: 12px !important;
        border-radius: 10px;
        transition: all 0.18s ease;
      }

      .navbar-nav .nav-link:hover,
      .navbar-nav .nav-link:focus {
        color: white !important;
        background: rgba(255,255,255,0.10);
      }

      .navbar-nav .nav-link.active,
      .navbar-nav .show > .nav-link {
        color: var(--ebi-green-dark) !important;
        background: white !important;
        box-shadow: 0 4px 14px rgba(0,0,0,0.10);
      }

      .bslib-sidebar,
      .bslib-sidebar .sidebar,
      .bslib-sidebar .sidebar-content,
      .bslib-sidebar .accordion,
      .bslib-sidebar .card,
      .bslib-sidebar > div {
        background: linear-gradient(180deg, #0B8F66 0%, #087556 100%) !important;
        border-color: rgba(255,255,255,0.12) !important;
        color: white !important;
      }

      .bslib-sidebar {
        border: 1px solid rgba(255,255,255,0.12) !important;
        box-shadow: var(--shadow-card);
      }

      .bslib-sidebar h4,
      .bslib-sidebar .form-label,
      .bslib-sidebar .control-label,
      .bslib-sidebar .subtle-note,
      .bslib-sidebar .shiny-text-output,
      .bslib-sidebar pre,
      .bslib-sidebar .form-check-label,
      .bslib-sidebar label,
      .bslib-sidebar .irs-min,
      .bslib-sidebar .irs-max,
      .bslib-sidebar .irs-from,
      .bslib-sidebar .irs-to,
      .bslib-sidebar .irs-single {
        color: white !important;
      }

      .bslib-sidebar hr {
        border-top: 1px solid rgba(255,255,255,0.18);
      }

      .bslib-sidebar .chip {
        background: rgba(255,255,255,0.14) !important;
        border: 1px solid rgba(255,255,255,0.18) !important;
        color: white !important;
        box-shadow: none;
      }

      .bslib-sidebar .chip b {
        color: #DFF7ED !important;
      }

      .bslib-sidebar .form-control,
      .bslib-sidebar .form-select {
        background: #FFFFFF !important;
        color: #0F172A !important;
        border: 1px solid rgba(255,255,255,0.18) !important;
        box-shadow: none !important;
      }

      .bslib-sidebar .form-check-input {
        background-color: #FFFFFF !important;
        border: 1px solid rgba(255,255,255,0.35) !important;
      }

      .bslib-sidebar .form-check-input:checked {
        background-color: #DFF7ED !important;
        border-color: #DFF7ED !important;
      }

      .bslib-sidebar .shiny-text-output,
      .bslib-sidebar pre {
        background: rgba(255,255,255,0.10) !important;
        border: 1px solid rgba(255,255,255,0.16) !important;
        border-radius: 14px;
        color: white !important;
      }

      .sidebar-scroll {
        position: sticky;
        top: 84px;
        height: calc(100vh - 104px);
        overflow-y: auto;
        padding-right: 10px;
        background: transparent !important;
      }

      .sidebar-scroll::-webkit-scrollbar { width: 8px; }
      .sidebar-scroll::-webkit-scrollbar-thumb {
        background: rgba(255,255,255,0.25);
        border-radius: 999px;
      }

      .main-fixed {
        height: calc(100vh - 90px);
        overflow: hidden;
        display: flex;
        align-items: stretch;
        justify-content: center;
        padding: 14px;
      }

      .page-wrap {
        height: calc(100vh - 90px);
        overflow-y: auto;
        padding: 18px;
      }

      .card {
        border: 1px solid rgba(15,23,42,0.06) !important;
        background: var(--panel) !important;
        box-shadow: var(--shadow-soft);
      }

      .card-header {
        background: transparent !important;
        border-bottom: 1px solid rgba(15,23,42,0.06) !important;
        font-weight: 800;
        color: var(--ink);
        font-size: 1rem;
        padding-top: 16px !important;
        padding-bottom: 14px !important;
      }

      .card-fullheight {
        height: 100%;
        width: 100%;
      }

      .chiprow {
        display: flex;
        gap: 12px;
        flex-wrap: wrap;
        margin-bottom: 14px;
      }

      .chip {
        padding: 10px 14px;
        border-radius: 999px;
        border: 1px solid rgba(11,143,102,0.12);
        background: linear-gradient(180deg, rgba(255,255,255,0.98), rgba(243,248,247,0.96));
        font-weight: 700;
        font-size: 0.92rem;
        color: var(--ink);
        box-shadow: 0 6px 18px rgba(15,23,42,0.05);
      }

      .chip b {
        font-weight: 900;
        color: var(--ebi-green);
      }

      .kpi-grid {
        display: grid;
        grid-template-columns: repeat(3, minmax(220px, 1fr));
        gap: 16px;
        margin-bottom: 18px;
      }

      .kpi-card {
        background: white;
        border: 1px solid rgba(11,143,102,0.10);
        border-radius: 18px;
        padding: 20px;
        box-shadow: 0 10px 24px rgba(15,23,42,0.05);
      }

      .kpi-label {
        color: var(--muted);
        font-weight: 700;
        font-size: 0.92rem;
        margin-bottom: 8px;
      }

      .kpi-value {
        color: var(--ebi-green-dark);
        font-size: 2rem;
        font-weight: 800;
        line-height: 1.1;
      }

      h3 {
        font-weight: 800;
        margin-bottom: 8px;
        color: var(--ink);
      }

      h4 {
        font-size: 1rem;
        font-weight: 800;
        color: var(--ink);
        margin-top: 6px;
        margin-bottom: 10px;
      }

      .subtle-note {
        color: var(--muted);
        font-size: 0.95rem;
        font-weight: 500;
      }

      hr {
        margin: 14px 0;
        border-top: 1px solid rgba(15,23,42,0.08);
      }

      .form-label, .control-label {
        font-weight: 700;
        color: var(--ink);
        margin-bottom: 6px;
      }

      .form-control, .form-select {
        border: 1px solid rgba(15,23,42,0.10) !important;
        box-shadow: none !important;
        background: rgba(255,255,255,0.95) !important;
      }

      .form-control:focus, .form-select:focus {
        border-color: var(--ebi-teal) !important;
        box-shadow: 0 0 0 0.2rem var(--ring) !important;
      }

      .form-check {
        margin-bottom: 8px;
      }

      table.dataTable {
        border-collapse: separate !important;
        border-spacing: 0 !important;
      }

      .dataTables_wrapper .dataTables_filter input,
      .dataTables_wrapper .dataTables_length select {
        border: 1px solid rgba(15,23,42,0.10) !important;
        border-radius: 10px !important;
        background: white !important;
        padding: 6px 10px !important;
      }

      table.dataTable thead th {
        background: #F3F8F7 !important;
        color: var(--ink) !important;
        border-bottom: 1px solid #DCE7E5 !important;
        font-weight: 800 !important;
      }

      table.dataTable.stripe tbody tr.odd,
      table.dataTable.display tbody tr.odd {
        background-color: rgba(248,252,251,0.75) !important;
      }

      table.dataTable tbody td {
        border-bottom: 1px solid rgba(220,231,229,0.8) !important;
      }

      pre, .shiny-text-output {
        background: rgba(255,255,255,0.95);
        border: 1px solid rgba(15,23,42,0.08);
        border-radius: 14px;
        padding: 12px 14px;
        color: var(--ink);
      }
    "))
  ),
  
  
  # TAB 1: CLUSTERS
 
  tabPanel(
    "Clusters",
    layout_sidebar(
      sidebar = sidebar(
        width = 420,
        div(
          class = "sidebar-scroll",
          div(
            class = "chiprow",
            div(class = "chip", HTML(paste0("Cells <b>", N_CELLS, "</b>"))),
            div(class = "chip", HTML(paste0("Clusters <b>", N_CLUST, "</b>")))
          ),
          
          tags$div(class = "subtle-note", "Interactive cluster embeddings for tSNE, UMAP, and PCA."),
          br(),
          
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
          checkboxInput("swap_xy", "Swap X / Y", value = FALSE),
          checkboxInput("flip_x", "Flip X axis", value = FALSE),
          checkboxInput("flip_y", "Flip Y axis", value = FALSE),
          
          hr(),
          h4("Colour and style"),
          sliderInput("pt_size", "Point size", min = 0.5, max = 6, value = 5, step = 0.1),
          sliderInput("pt_alpha", "Point alpha", min = 0.1, max = 1.0, value = 0.88, step = 0.05),
          checkboxInput("show_labels", "Show cluster labels", value = TRUE),
          checkboxInput("show_legend", "Show legend", value = TRUE),
          
          hr(),
          h4("Run summary"),
          verbatimTextOutput("run_summary")
        )
      ),
      
      div(
        class = "main-fixed",
        card(
          class = "card-fullheight",
          card_header(
            div(
              style = "display:flex; justify-content:space-between; align-items:center;",
              tags$span("Cluster map")
            )
          ),
          plotOutput("cluster_plot", height = "calc(100vh - 170px)")
        )
      )
    )
  ),
  
 
  # TAB 2: GENE VIEW

  tabPanel(
    "Gene View",
    layout_sidebar(
      sidebar = sidebar(
        width = 420,
        div(
          class = "sidebar-scroll",
          div(
            class = "chiprow",
            div(class = "chip", "Fallback mode"),
            div(class = "chip", "Gene view disabled")
          ),
          
          h4("Status"),
          tags$div(
            class = "subtle-note",
            "This fallback build runs without Seurat and without expression matrix support."
          ),
          br(),
          tags$div(
            class = "subtle-note",
            "Gene-level plots and gene summary tables are disabled in this mode."
          )
        )
      ),
      div(
        class = "main-fixed",
        div(
          style = "width:100%; height:100%; overflow:hidden;",
          card(
            class = "card-fullheight",
            card_header(
              div(
                style = "display:flex; justify-content:space-between; align-items:center;",
                tags$span("Gene View")
              )
            ),
            div(
              style = "padding: 32px; font-size: 1rem; line-height: 1.8;",
              tags$h4("Fallback mode active"),
              tags$p("This app version preserves the same UI design but removes Seurat and Plotly dependencies."),
              tags$p("The Clusters tab is fully functional using the precomputed embedding files."),
              tags$p("The Gene View tab is shown only as a placeholder because expression data is not available in this mode.")
            )
          )
        )
      )
    )
  ),
  
 
  # TAB 3: SUMMARY
 
  tabPanel(
    "Summary",
    div(
      class = "page-wrap",
      h3("Dataset Summary"),
      tags$div(
        class = "subtle-note",
        "Overview of the loaded single-cell dataset and precomputed dimensionality reductions."
      ),
      br(),
      
      div(
        class = "kpi-grid",
        div(
          class = "kpi-card",
          div(class = "kpi-label", "Total cells"),
          div(class = "kpi-value", N_CELLS)
        ),
        div(
          class = "kpi-card",
          div(class = "kpi-label", "Clusters"),
          div(class = "kpi-value", N_CLUST)
        ),
        div(
          class = "kpi-card",
          div(class = "kpi-label", "Mode"),
          div(class = "kpi-value", "Fallback")
        )
      ),
      
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Project summary"),
          tags$div(
            style = "padding: 6px 2px 4px 2px; line-height:1.8;",
            tags$p("This app presents a fast interactive explorer for the GSE75140 precomputed embeddings."),
            tags$p("The Clusters tab allows inspection of precomputed tSNE, UMAP, and PCA embeddings with alignment controls."),
            tags$p("This fallback version is designed to run even when Seurat and Plotly are unavailable."),
            tags$p("Gene-level expression analysis is disabled in this build because no expression matrix is loaded.")
          )
        ),
        card(
          card_header("Loaded components"),
          tags$div(
            style = "padding: 6px 2px 4px 2px; line-height:1.8;",
            tags$p(HTML(paste0("<b>tSNE embedding:</b> ", TSNE_RDS))),
            tags$p(HTML(paste0("<b>UMAP embedding:</b> ", UMAP_RDS))),
            tags$p(HTML(paste0("<b>PCA embedding:</b> ", PCA_RDS))),
            tags$p(HTML(paste0("<b>Species metadata available:</b> ", ifelse(species_exists, "Yes", "No"))))
          )
        )
      ),
      
      br(),
      
      card(
        card_header("Cluster identities"),
        DTOutput("cluster_table")
      )
    )
  )
)


# Server

server <- function(input, output, session) {
  
  output$run_summary <- renderText({
    paste0(
      "Mode: fallback / no Seurat / no plotly\n",
      "Cells loaded: ", N_CELLS, "\n",
      "Clusters loaded: ", N_CLUST, "\n",
      "Embeddings available: tSNE, UMAP, PCA\n",
      "Gene View: placeholder only"
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
    df <- safe_species(df)
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
      scale_color_manual(values = CLUST_COLORS, drop = FALSE) +
      labs(
        x = input$x_axis,
        y = input$y_axis,
        color = "Cluster"
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", colour = "#0F172A"),
        legend.position = if (isTRUE(input$show_legend)) "right" else "none",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "#0F172A")
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
          size = 4.8,
          fontface = "bold"
        )
    }
    
    p
  })
  
  output$cluster_table <- renderDT({
    cl <- data.frame(
      cluster = cluster_levels,
      stringsAsFactors = FALSE
    )
    
    datatable(
      cl,
      options = list(pageLength = 12, dom = "tip"),
      rownames = FALSE
    )
  })
}

shinyApp(ui, server)
