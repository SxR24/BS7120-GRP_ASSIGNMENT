library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(plotly)
library(dplyr)
library(Seurat)

# Precomputed files

SEURAT_RDS <- "gse75140_seurat_ready.rds"
TSNE_RDS   <- "emb_tsne.rds"
UMAP_RDS   <- "emb_umap.rds"
PCA_RDS    <- "emb_pca.rds"



make_palette <- function(levels_vec) {
  cols <- grDevices::hcl.colors(length(levels_vec), palette = "Dark 3")
  names(cols) <- levels_vec
  cols
}

get_expr_mat <- function(seu_obj) {
  a <- seu_obj[["RNA"]]
  lyr <- Layers(a)
  if ("data" %in% lyr) {
    return(GetAssayData(seu_obj, layer = "data"))
  }
  GetAssayData(seu_obj, layer = "counts")
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


seurat_obj <- readRDS(SEURAT_RDS)
emb_tsne   <- readRDS(TSNE_RDS)
emb_umap   <- readRDS(UMAP_RDS)
emb_pca    <- readRDS(PCA_RDS)

expr_mat <- get_expr_mat(seurat_obj)

ALL_GENES <- rownames(seurat_obj)
cluster_levels <- sort(unique(as.character(Idents(seurat_obj))))
CLUST_COLORS <- make_palette(cluster_levels)

# Basic run summary values
N_CELLS <- ncol(seurat_obj)
N_GENES <- nrow(seurat_obj)
N_CLUST <- length(unique(Idents(seurat_obj)))

# Theme

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

# UI

ui <- navbarPage(
  title = tagList(
    tags$span(style = "font-weight:700; letter-spacing:0.2px;", "GSE75140 Explorer"),
    tags$span(style = "opacity:0.7; font-weight:600; margin-left:8px; font-size:0.92rem;", "• fast mode")
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

      .navbar .navbar-brand span { color: var(--nav-muted) !important; }

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

      .navbar-toggler {
        border-color: rgba(15,23,42,0.20) !important;
      }

      .navbar-toggler-icon {
        filter: invert(0) !important;
        opacity: 0.85;
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
      .plotly { height: 100% !important; }

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
            div(class = "chip", HTML(paste0("Genes: <b>", N_GENES, "</b>"))),
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
          h4("Colour & style"),
          sliderInput("pt_size", "Point size", min = 1, max = 12, value = 3, step = 1),
          sliderInput("pt_alpha", "Point alpha", min = 0.1, max = 1.0, value = 0.9, step = 0.05),
          checkboxInput("show_labels", "Show labels (c1..)", value = TRUE),
          checkboxInput("show_legend", "Show legend", value = TRUE),
          
          hr(),
          h4("Gene overlay"),
          checkboxInput("use_gene_overlay", "Overlay gene expression", value = FALSE),
          selectizeInput(
            "gene_overlay", "Gene",
            choices = NULL,
            options = list(placeholder = "Type gene...", maxOptions = 5000)
          ),
          checkboxInput("log_expr", "Log-scale (log1p)", value = TRUE),
          
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
          plotlyOutput("cluster_plot", height = "calc(100vh - 170px)")
        )
      )
    )
  ),
  
  tabPanel(
    "Gene View",
    layout_sidebar(
      sidebar = sidebar(
        width = 430,
        div(
          class = "sidebar-scroll",
          h4("Gene query"),
          selectizeInput(
            "gene", "Gene",
            choices = NULL,
            options = list(placeholder = "Type gene...", maxOptions = 5000)
          ),
          selectInput(
            "group_by", "Group by",
            choices = c("Cluster" = "cluster", "Species" = "species"),
            selected = "cluster"
          ),
          checkboxInput("violin_log", "Log-scale violin (log1p)", value = TRUE)
        )
      ),
      div(
        class = "main-fixed",
        div(
          style = "width:100%; height:100%; overflow:hidden;",
          layout_columns(
            col_widths = c(6, 6),
            card(
              class = "card-fullheight",
              card_header("Violin"),
              plotOutput("violin_plot", height = "calc(100vh - 320px)")
            ),
            card(
              class = "card-fullheight",
              card_header("Summary"),
              DTOutput("gene_summary")
            )
          ),
          card(
            class = "card-fullheight",
            card_header("Top expressing cells (top 200)"),
            DTOutput("gene_cells")
          )
        )
      )
    )
  )
)


# Server

server <- function(input, output, session) {
  
  # One-time initialization for gene dropdowns
  session$onFlushed(function() {
    default_gene <- if ("FOXG1" %in% ALL_GENES) "FOXG1" else ALL_GENES[1]
    
    updateSelectizeInput(
      session, "gene_overlay",
      choices = ALL_GENES,
      selected = default_gene,
      server = TRUE
    )
    
    updateSelectizeInput(
      session, "gene",
      choices = ALL_GENES,
      selected = default_gene,
      server = TRUE
    )
  }, once = TRUE)
  
  output$run_summary <- renderText({
    paste0(
      "Mode: precomputed / fast load\n",
      "Cells loaded: ", N_CELLS, "\n",
      "Genes loaded: ", N_GENES, "\n",
      "Clusters loaded: ", N_CLUST, "\n",
      "Active assay: RNA\n",
      "Embeddings available: tSNE, UMAP, PCA"
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
    num_cols <- names(df)[sapply(df, is.numeric)]
    
    preferred <- num_cols[grepl("^(tSNE|UMAP|PC)_", num_cols)]
    if (length(preferred) >= 2) {
      choices <- preferred
      x_default <- preferred[1]
      y_default <- preferred[2]
    } else {
      choices <- num_cols
      x_default <- num_cols[1]
      y_default <- num_cols[2]
    }
    
    tagList(
      selectInput("x_axis", "X axis", choices = choices, selected = x_default),
      selectInput("y_axis", "Y axis", choices = choices, selected = y_default)
    )
  })
  
  plot_df <- reactive({
    req(input$x_axis, input$y_axis)
    
    df <- embedding_df()
    df <- apply_alignment(
      df,
      input$x_axis,
      input$y_axis,
      input$swap_xy,
      input$flip_x,
      input$flip_y
    )
    
    df$cluster <- factor(df$cluster, levels = names(CLUST_COLORS))
    df
  })
  
  output$cluster_plot <- renderPlotly({
    df <- plot_df()
    
    axis_common <- list(
      showgrid  = TRUE,
      zeroline  = FALSE,
      showline  = TRUE,
      linecolor = "rgba(0,0,0,0.25)",
      ticks     = "outside",
      tickcolor = "rgba(0,0,0,0.25)",
      ticklen   = 5
    )
    
    if (isTRUE(input$use_gene_overlay)) {
      gene <- input$gene_overlay
      req(gene, gene %in% ALL_GENES)
      
      v <- as.numeric(expr_mat[gene, df$cell_id])
      if (isTRUE(input$log_expr)) v <- log1p(v)
      df$expr <- v
      
      p <- plot_ly(
        data = df,
        x = ~X_plot, y = ~Y_plot,
        type = "scattergl",
        mode = "markers",
        marker = list(size = input$pt_size, opacity = input$pt_alpha),
        color = ~expr,
        text = ~paste0(
          "cell: ", cell_id, "<br>",
          "species: ", species, "<br>",
          "cluster: ", as.character(cluster), "<br>",
          gene, ": ", round(expr, 3)
        ),
        hoverinfo = "text"
      ) %>%
        layout(
          plot_bgcolor = "rgba(255,255,255,0)",
          paper_bgcolor = "rgba(255,255,255,0)",
          xaxis = c(list(title = input$x_axis), axis_common),
          yaxis = c(list(title = input$y_axis), axis_common),
          showlegend = FALSE
        )
      
      return(p)
    }
    
    p <- plot_ly(
      data = df,
      x = ~X_plot, y = ~Y_plot,
      type = "scattergl",
      mode = "markers",
      marker = list(size = input$pt_size, opacity = input$pt_alpha),
      color = ~cluster,
      colors = CLUST_COLORS,
      text = ~paste0(
        "cell: ", cell_id, "<br>",
        "species: ", species, "<br>",
        "cluster: ", as.character(cluster)
      ),
      hoverinfo = "text"
    ) %>%
      layout(
        plot_bgcolor = "rgba(255,255,255,0)",
        paper_bgcolor = "rgba(255,255,255,0)",
        xaxis = c(list(title = input$x_axis), axis_common),
        yaxis = c(list(title = input$y_axis), axis_common),
        showlegend = isTRUE(input$show_legend),
        legend = list(title = list(text = "Cluster"))
      )
    
    if (isTRUE(input$show_labels)) {
      cents <- cluster_centroids(df, "cluster")
      cents$label <- as.character(cents$cluster)
      
      p <- p %>%
        add_text(
          data = cents,
          x = ~Xc, y = ~Yc,
          text = ~label,
          textfont = list(color = "black", size = 13),
          inherit = FALSE,
          showlegend = FALSE
        )
    }
    
    p
  })
  
  gene_df_for_plot <- reactive({
    gene <- input$gene
    req(gene, gene %in% ALL_GENES)
    
    cells <- colnames(expr_mat)
    v <- as.numeric(expr_mat[gene, cells])
    
    if (isTRUE(input$violin_log)) {
      v <- log1p(v)
    }
    
    data.frame(
      cell_id = cells,
      expr = v,
      cluster = factor(as.character(Idents(seurat_obj)[cells]), levels = names(CLUST_COLORS)),
      species = seurat_obj@meta.data[cells, "species", drop = TRUE],
      stringsAsFactors = FALSE
    )
  })
  
  output$violin_plot <- renderPlot({
    df <- gene_df_for_plot()
    gene <- input$gene
    
    if (input$group_by == "cluster") {
      ggplot(df, aes(x = cluster, y = expr, fill = cluster)) +
        geom_violin(trim = TRUE) +
        geom_jitter(width = 0.2, size = 0.5, alpha = 0.25) +
        scale_fill_manual(values = CLUST_COLORS) +
        labs(
          title = paste0("Expression: ", gene),
          x = "Cluster",
          y = ifelse(input$violin_log, "log1p(expr)", "expr")
        ) +
        theme_classic(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none"
        )
    } else {
      ggplot(df, aes(x = species, y = expr)) +
        geom_violin(trim = TRUE, fill = "#93c5fd") +
        geom_jitter(width = 0.2, size = 0.5, alpha = 0.25) +
        labs(
          title = paste0("Expression: ", gene),
          x = "Species",
          y = ifelse(input$violin_log, "log1p(expr)", "expr")
        ) +
        theme_classic(base_size = 12) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
    }
  })
  
  output$gene_summary <- renderDT({
    df <- gene_df_for_plot()
    
    if (input$group_by == "cluster") {
      out <- df %>%
        group_by(cluster) %>%
        summarise(
          n_cells = n(),
          pct_expr_gt0 = mean(expr > 0) * 100,
          mean_expr = mean(expr),
          median_expr = median(expr),
          .groups = "drop"
        ) %>%
        arrange(desc(mean_expr))
    } else {
      out <- df %>%
        group_by(species) %>%
        summarise(
          n_cells = n(),
          pct_expr_gt0 = mean(expr > 0) * 100,
          mean_expr = mean(expr),
          median_expr = median(expr),
          .groups = "drop"
        ) %>%
        arrange(desc(mean_expr))
    }
    
    datatable(
      out,
      options = list(pageLength = 12, scrollX = TRUE)
    )
  })
  
  output$gene_cells <- renderDT({
    df <- gene_df_for_plot()
    
    top <- df %>%
      arrange(desc(expr)) %>%
      head(200)
    
    datatable(
      top,
      options = list(pageLength = 12, scrollX = TRUE)
    )
  })
}

shinyApp(ui, server)
