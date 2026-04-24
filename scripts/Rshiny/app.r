library(shiny)

# Plot file

original_plots <- list(
  "Figure S1A" = c("Figure_S1A.rds"),
  "Figure 1B Heatmap A" = c("Figure_1B_heatmapa.rds"),
  "Figure 1B Heatmap B" = c("Figure_1B_heatmapb.rds"),
  "Figure 1C Heatmap" = c("Figure_1C_heatmap_plot.rds"),
  "Figure 1C Heatmap Part 2" = c("Figure_1C_part2_heatmap_plot.rds"),
  "Figure S1B Heatmap" = c("Figure_S1B_heatmap_plot.rds"),
  "Figure S1C Heatmap" = c("Figure_S1C_heatmap_plot.rds"),
  "Original Figure 2A" = c("OldpipelineFig2A.rds"),
  "Original Figure 2B.1" = c("OldpipelineFig2B1.rds"),
  "Original Figure 2B.2" = c("OldpipelineFig2B2.rds"),
  "Original Figure 2B.3" = c("OldpipelineFig2B3.rds"),
  "Original Figure 2B.4" = c("OldpipelineFig2B4.rds"),
  "Figure 3D" = c("Fig3D_plot.rds"),
  "Original Figure 3E" = c("OldpipelineFig3E.rds"),
  "Figure S3H" = c("FigS3H.rds"),
  "Organoid Cluster Comparison" = c("Organoidclustercomparison.rds"),
  "PCA Overlap Plot" = c("PCA_overlap_plot.rds"),
  "Figure 1C Comparison" = c("Fig1Ccomparison.rds"),
  "Figure S1A Comparison" = c("FigS1Acomparison.rds"),
  "Figure S1A Comparison 2" = c("FigS1Acomparison2.rds", "Fig1SAcomparison2.rds")
)

alternative_plots <- list(
  "Alternative Figure 3D" = c("Alt_Fig3D.rds", "Alt_Fig3D - Copy.rds"),
  "Alternative Figure S3H" = c("Alt_FigS3H.rds", "Alt_FigS3H - Copy.rds"),
  "Protein Prefilter" = c("Altpipelineprotienprefilter.rds", "Altpipelineproteinprefilter.rds"),
  "Protein Figure S1A" = c("AltFigure_S1A_protein.rds", "AltFigure_S1A_protein - Copy.rds"),
  "Protein Figure S1B Heatmap" = c(
    "AltFigure_S1B_heatmap_plot_protein.rds",
    "AltFigure_S1B_heatmap_plot_protein - Copy.rds",
    "AltFigure_S1B_heatmap_plot_protien.rds",
    "AltFigure_S1B_heatmap_plot_protien - Copy.rds"
  ),
  "Protein Figure S1C Heatmap" = c(
    "AltFigure_S1c_heatmap_plot_protein.rds",
    "AltFigure_S1c_heatmap_plot_protein - Copy.rds",
    "AltFigure_S1c_heatmap_plot_protien.rds",
    "AltFigure_S1c_heatmap_plot_protien - Copy.rds"
  ),
  "Protein Figure 3E" = c("AltprotienFig3E.rds"),
  "Protein Cell-Type Heatmap" = c(
    "Altpipeline_protien_fetalcelltypeheatmap.rds",
    "Altpipeline_protein_fetalcelltypeheatmap.rds"
  ),
  "Protein Cluster Expression" = c(
    "Altpipeline_protien_fetalclusterexpression.rds",
    "Altpipeline_protein_fetalclusterexpression.rds"
  ),
  "Protein Cluster-Zone Correlation" = c(
    "Altpipeline_protien_fetalclusterzonecorrelation.rds",
    "Altpipeline_protein_fetalclusterzonecorrelation.rds"
  ),
  "Protein Proto-Cluster" = c(
    "Altpipeline_protien_fetalprotocluster.rds",
    "Altpipeline_protein_fetalprotocluster.rds"
  ),
  "Protein Subcluster" = c(
    "Altpipeline_protien_fetalsubcluster.rds",
    "Altpipeline_protein_fetalsubcluster.rds"
  ),
  "Protein Subcluster Expression" = c(
    "Altpipeline_protien_fetalsubclusterexpression.rds",
    "Altpipeline_protein_fetalsubclusterexpression.rds"
  ),
  "Protein Subpopulation Cluster" = c(
    "Altpipeline_protien_fetalsubpoplationcluster.rds",
    "Altpipeline_protein_fetalsubpoplationcluster.rds"
  ),
  "Protein Zone Heatmap" = c(
    "Altpipeline_protien_fetalzoneheatmap.rds",
    "Altpipeline_protein_fetalzoneheatmap.rds"
  ),
  "Non-Coding Prefilter" = c("Altpipelinenoncodingprefilter.rds"),
  "Non-Coding Figure S1A" = c("AltFigure_S1A_noncoding.rds", "AltFigure_S1A_noncoding - Copy.rds"),
  "Non-Coding Figure S1B Heatmap" = c("AltFigure_S1B_heatmap_plot_noncoding.rds", "AltFigure_S1B_heatmap_plot_noncoding - Copy.rds"),
  "Non-Coding Figure S1C Heatmap" = c("AltFigure_S1c_heatmap_plot_noncoding.rds", "AltFigure_S1c_heatmap_plot_noncoding - Copy.rds"),
  "Non-Coding Figure 3D" = c("AltnoncodingFig3D.rds", "AltnoncodingFig3D - Copy.rds"),
  "Non-Coding Figure 3E" = c("AltnoncodingFig3E.rds", "AltnoncodingFig3E - Copy.rds"),
  "Non-Coding Figure S3H" = c("AltnoncodingFigS3H.rds", "AltnoncodingFigS3H - Copy.rds"),
  "Non-Coding Cell-Type Heatmap" = c("Altpipeline_noncoding_fetalcelltypeheatmap.rds"),
  "Non-Coding Cluster Expression" = c("Altpipeline_noncoding_fetalclusterexpression.rds"),
  "Non-Coding Cluster-Zone Correlation" = c("Altpipeline_noncoding_fetalclusterzonecorrelation.rds"),
  "Non-Coding Proto-Cluster" = c("Altpipeline_noncoding_fetalprotocluster.rds"),
  "Non-Coding Subcluster" = c("Altpipeline_noncoding_fetalsubcluster.rds"),
  "Non-Coding Subcluster Expression" = c("Altpipeline_noncoding_fetalsubclusterexpression.rds"),
  "Non-Coding Subpopulation Cluster" = c("Altpipeline_noncoding_fetalsubpoplationcluster.rds"),
  "Non-Coding Zone Heatmap" = c("Altpipeline_noncoding_fetalzoneheatmap.rds")
)

# Annotation

plot_notes <- list(
  "Figure S1A" = "Write annotation here.",
  "Figure 1B Heatmap A" = "Write annotation here.",
  "Figure 1B Heatmap B" = "Write annotation here.",
  "Figure 1C Heatmap" = "Write annotation here.",
  "Figure 1C Heatmap Part 2" = "Write annotation here.",
  "Figure S1B Heatmap" = "Write annotation here.",
  "Figure S1C Heatmap" = "Write annotation here.",
  "Original Figure 2A" = "Write annotation here.",
  "Original Figure 2B.1" = "Write annotation here.",
  "Original Figure 2B.2" = "Write annotation here.",
  "Original Figure 2B.3" = "Write annotation here.",
  "Original Figure 2B.4" = "Write annotation here.",
  "Figure 3D" = "Write annotation here.",
  "Original Figure 3E" = "Write annotation here.",
  "Figure S3H" = "Write annotation here.",
  "Organoid Cluster Comparison" = "Write annotation here.",
  "PCA Overlap Plot" = "Write annotation here.",
  "Figure 1C Comparison" = "Write annotation here.",
  "Figure S1A Comparison" = "Write annotation here.",
  "Figure S1A Comparison 2" = "Write annotation here.",
  "Alternative Figure 3D" = "Write annotation here.",
  "Alternative Figure S3H" = "Write annotation here.",
  "Protein Prefilter" = "Write annotation here.",
  "Protein Figure S1A" = "Write annotation here.",
  "Protein Figure S1B Heatmap" = "Write annotation here.",
  "Protein Figure S1C Heatmap" = "Write annotation here.",
  "Protein Figure 3E" = "Write annotation here.",
  "Protein Cell-Type Heatmap" = "Write annotation here.",
  "Protein Cluster Expression" = "Write annotation here.",
  "Protein Cluster-Zone Correlation" = "Write annotation here.",
  "Protein Proto-Cluster" = "Write annotation here.",
  "Protein Subcluster" = "Write annotation here.",
  "Protein Subcluster Expression" = "Write annotation here.",
  "Protein Subpopulation Cluster" = "Write annotation here.",
  "Protein Zone Heatmap" = "Write annotation here.",
  "Non-Coding Prefilter" = "Write annotation here.",
  "Non-Coding Figure S1A" = "Write annotation here.",
  "Non-Coding Figure S1B Heatmap" = "Write annotation here.",
  "Non-Coding Figure S1C Heatmap" = "Write annotation here.",
  "Non-Coding Figure 3D" = "Write annotation here.",
  "Non-Coding Figure 3E" = "Write annotation here.",
  "Non-Coding Figure S3H" = "Write annotation here.",
  "Non-Coding Cell-Type Heatmap" = "Write annotation here.",
  "Non-Coding Cluster Expression" = "Write annotation here.",
  "Non-Coding Cluster-Zone Correlation" = "Write annotation here.",
  "Non-Coding Proto-Cluster" = "Write annotation here.",
  "Non-Coding Subcluster" = "Write annotation here.",
  "Non-Coding Subcluster Expression" = "Write annotation here.",
  "Non-Coding Subpopulation Cluster" = "Write annotation here.",
  "Non-Coding Zone Heatmap" = "Write annotation here."
)

# Axis Inoformation

interactive_axis_info <- list(
  "pca" = data.frame(
    Axis = c("X Axis", "Y Axis"),
    Description = c("PCA or overlap axis 1", "PCA or overlap axis 2"),
    stringsAsFactors = FALSE
  ),
  "org" = data.frame(
    Axis = c("X Axis", "Y Axis"),
    Description = c("Organoid comparison axis 1", "Organoid comparison axis 2"),
    stringsAsFactors = FALSE
  ),
  "orig3d" = data.frame(
    Axis = c("X Axis", "Y Axis"),
    Description = c("Original Figure 3D tSNE 1", "Original Figure 3D tSNE 2"),
    stringsAsFactors = FALSE
  ),
  "alt3d" = data.frame(
    Axis = c("X Axis", "Y Axis"),
    Description = c("Alternative Figure 3D axis 1", "Alternative Figure 3D axis 2"),
    stringsAsFactors = FALSE
  ),
  "nc3d" = data.frame(
    Axis = c("X Axis", "Y Axis"),
    Description = c("Non coding Figure 3D axis 1", "Non coding Figure 3D axis 2"),
    stringsAsFactors = FALSE
  )
)

# Conclusion

conclusion_text <- ""

# Data files

cell_file <- if (file.exists("cell_assignments.txt")) "cell_assignments.txt" else "cell_assignments"
cell_data <- read.table(cell_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

gene_matrix_file <- if (file.exists("normalized_expression.rds")) {
  "normalized_expression.rds"
} else if (file.exists("/mnt/data/normalized_expression.rds")) {
  "/mnt/data/normalized_expression.rds"
} else {
  stop("Gene matrix file was not found.")
}

# User Interface UI

ui_base <- fluidPage(
  tags$head(
    tags$style(HTML("
      html, body { height: 100%; margin: 0; background: #f5f7fb; font-family: 'Segoe UI', Arial, sans-serif; color: #172033; }
      .container-fluid { padding-left: 0; padding-right: 0; }
      .topbar { min-height: 68px; background: linear-gradient(90deg, #02172f 0%, #041a3c 100%); display: flex; align-items: center; justify-content: space-between; padding: 0 22px; box-shadow: 0 2px 8px rgba(0,0,0,0.12); }
      .topbar-left { display: flex; align-items: center; gap: 20px; flex-wrap: wrap; }
      .brand { color: white; font-size: 18px; font-weight: 700; display: flex; align-items: center; gap: 12px; }
      .brand-icon { width: 26px; height: 26px; border-radius: 50%; border: 2px solid #1e90ff; display: inline-flex; align-items: center; justify-content: center; color: #1e90ff; font-weight: 700; }
      .top-nav { display: flex; align-items: center; gap: 10px; flex-wrap: wrap; }
      .nav-btn { background: transparent; color: #e7edf8; border: none; border-radius: 12px; padding: 12px 18px; font-size: 15px; font-weight: 600; cursor: pointer; }
      .nav-btn.active-tab { background: #153b79; color: white; }

      .app-shell { display: flex; height: calc(100vh - 68px); }
      .sidebar { width: 340px; min-width: 340px; background: #f8fafc; border-right: 1px solid #e2e8f0; padding: 18px 16px; box-sizing: border-box; }
      .sidebar-card { background: #ffffff; border-radius: 16px; box-shadow: 0 3px 14px rgba(15, 23, 42, 0.07); border: 1px solid #e5eaf2; padding: 16px; height: calc(100vh - 104px); display: flex; flex-direction: column; }
      .side-section-title { font-size: 12px; font-weight: 800; letter-spacing: 0.6px; color: #25406e; margin-bottom: 12px; text-transform: uppercase; }
      .plot-select-wrap { margin-top: 8px; }
      .plot-select-wrap .form-control { border-radius: 12px; border: 1px solid #d8e0ee; height: 46px; box-shadow: none; }
      .plot-list-wrap { margin-top: 14px; flex: 1 1 auto; min-height: 0; overflow-y: auto; border-top: 1px solid #eef2f7; padding-top: 10px; }
      .plot-list-button { width: 100%; text-align: left; border: none; background: transparent; border-radius: 14px; padding: 14px 14px; margin-bottom: 8px; color: #24344d; font-size: 15px; font-weight: 600; cursor: pointer; }
      .plot-list-button:hover { background: #eef4ff; }
      .plot-list-button.active-plot { background: #dbe6fb; color: #1c49b8; font-weight: 700; }

      .main-content { flex: 1; padding: 20px 22px; overflow: auto; }
      .content-header { display: flex; align-items: flex-start; justify-content: space-between; gap: 20px; margin-bottom: 16px; }
      .content-title { font-size: 28px; font-weight: 800; color: #172033; margin: 0 0 6px 0; }
      .breadcrumb-line { font-size: 15px; color: #5a6c8f; }
      .breadcrumb-line .crumb-accent { color: #2352d7; font-weight: 600; }

      .plot-panel { background: white; border-radius: 16px; border: 1px solid #e5eaf2; box-shadow: 0 3px 14px rgba(15, 23, 42, 0.07); padding: 8px; margin-bottom: 16px; }
      .annotation-panel { background: white; border-radius: 16px; border: 1px solid #e5eaf2; box-shadow: 0 3px 14px rgba(15, 23, 42, 0.07); padding: 22px; position: relative; margin-bottom: 16px; }
      .annotation-panel:before { content: ''; position: absolute; left: 0; top: 28px; bottom: 28px; width: 4px; background: #2563eb; border-radius: 3px; }
      .annotation-title { font-size: 16px; font-weight: 800; color: #172033; margin-bottom: 14px; padding-left: 18px; }
      .annotation-text { color: #4a5a75; font-size: 15px; line-height: 1.8; padding-left: 18px; white-space: pre-wrap; }

      .interactive-shell { padding: 20px 22px; }
      .interactive-card { background: white; border-radius: 16px; border: 1px solid #e5eaf2; box-shadow: 0 3px 14px rgba(15, 23, 42, 0.07); padding: 18px; margin-bottom: 16px; }
      .control-label { font-weight: 700; color: #334155; }
      .form-control { border-radius: 12px; border: 1px solid #d8e0ee; box-shadow: none; }
      .table { background: white; }

      .conclusion-box { background: white; border-radius: 16px; border: 1px solid #e5eaf2; box-shadow: 0 3px 14px rgba(15, 23, 42, 0.07); padding: 30px; min-height: 300px; }

      .axis-card, .gene-card { background: white; border-radius: 16px; border: 1px solid #e5eaf2; box-shadow: 0 3px 14px rgba(15, 23, 42, 0.07); padding: 18px; margin-top: 16px; }

      .axis-table-wrap { width: 100%; overflow-x: auto; }
      .axis-table { width: 100%; border-collapse: separate; border-spacing: 0; font-size: 15px; color: #24344d; background: #ffffff; border: 1px solid #e2e8f0; border-radius: 12px; overflow: hidden; }
      .axis-table thead th { background: #eff4fb; color: #10233f; font-weight: 700; text-align: left; padding: 14px 16px; border-bottom: 1px solid #dbe5f0; }
      .axis-table tbody td { padding: 14px 16px; border-bottom: 1px solid #edf2f7; vertical-align: top; }
      .axis-table tbody tr:last-child td { border-bottom: none; }
      .axis-table tbody tr:nth-child(even) { background: #fafcff; }
      .axis-col { width: 28%; font-weight: 700; color: #1c49b8; white-space: nowrap; }
      .desc-col { width: 72%; color: #334155; }

      .gene-layout { padding: 20px 22px; }
      .gene-page-grid { display: grid; grid-template-columns: 340px 1fr; gap: 18px; align-items: start; }
      .gene-control-card { background: white; border-radius: 16px; border: 1px solid #e5eaf2; box-shadow: 0 3px 14px rgba(15, 23, 42, 0.07); padding: 18px; }
      .gene-result-card { background: white; border-radius: 16px; border: 1px solid #e5eaf2; box-shadow: 0 3px 14px rgba(15, 23, 42, 0.07); padding: 18px; margin-bottom: 16px; }
      .small-muted { color: #64748b; font-size: 14px; line-height: 1.6; }
    "))
  ),
  
  div(
    class = "topbar",
    div(
      class = "topbar-left",
      div(class = "brand", span(class = "brand-icon", "◉"), "Group-A Dashboard Analysis of Human Cerebral Organoids"),
      div(
        class = "top-nav",
        actionButton("nav_original", "Original Pipeline", class = "nav-btn"),
        actionButton("nav_alternative", "Alternative Pipeline", class = "nav-btn"),
        actionButton("nav_interactive", "Interactive Panel", class = "nav-btn"),
        actionButton("nav_gene", "Gene Query", class = "nav-btn"),
        actionButton("nav_conclusion", "Conclusion", class = "nav-btn")
      )
    ),
    div()
  ),
  
  uiOutput("page_ui")
)

# Server 

server <- function(input, output, session) {
  
  plot_cache <- reactiveVal(new.env(parent = emptyenv()))
  gene_data_cache <- reactiveVal(NULL)
  current_page <- reactiveVal("original")
  
  observeEvent(input$nav_original, current_page("original"))
  observeEvent(input$nav_alternative, current_page("alternative"))
  observeEvent(input$nav_interactive, current_page("interactive"))
  observeEvent(input$nav_gene, current_page("gene"))
  observeEvent(input$nav_conclusion, current_page("conclusion"))
  
  output$orig_plot_list_ui <- renderUI({
    req(input$orig_plot_choice)
    tagList(lapply(names(original_plots), function(p) {
      btn_id <- paste0("orig_plot_btn_", gsub("[^A-Za-z0-9]", "_", p))
      cls <- if (identical(p, input$orig_plot_choice)) "plot-list-button active-plot" else "plot-list-button"
      actionButton(btn_id, label = p, class = cls)
    }))
  })
  
  output$alt_plot_list_ui <- renderUI({
    req(input$alt_plot_choice)
    tagList(lapply(names(alternative_plots), function(p) {
      btn_id <- paste0("alt_plot_btn_", gsub("[^A-Za-z0-9]", "_", p))
      cls <- if (identical(p, input$alt_plot_choice)) "plot-list-button active-plot" else "plot-list-button"
      actionButton(btn_id, label = p, class = cls)
    }))
  })
  
  observe({
    lapply(names(original_plots), function(p) {
      local({
        plot_name <- p
        btn_id <- paste0("orig_plot_btn_", gsub("[^A-Za-z0-9]", "_", plot_name))
        observeEvent(input[[btn_id]], {
          updateSelectInput(session, "orig_plot_choice", selected = plot_name)
        }, ignoreInit = TRUE)
      })
    })
  })
  
  observe({
    lapply(names(alternative_plots), function(p) {
      local({
        plot_name <- p
        btn_id <- paste0("alt_plot_btn_", gsub("[^A-Za-z0-9]", "_", plot_name))
        observeEvent(input[[btn_id]], {
          updateSelectInput(session, "alt_plot_choice", selected = plot_name)
        }, ignoreInit = TRUE)
      })
    })
  })
  
  output$page_ui <- renderUI({
    if (current_page() == "original") {
      div(
        class = "app-shell",
        div(
          class = "sidebar",
          div(
            class = "sidebar-card",
            div(class = "side-section-title", "Original Pipeline"),
            div(
              class = "plot-select-wrap",
              selectInput(
                "orig_plot_choice",
                NULL,
                choices = names(original_plots),
                selected = if (is.null(input$orig_plot_choice)) names(original_plots)[1] else input$orig_plot_choice,
                width = "100%"
              )
            ),
            div(class = "plot-list-wrap", uiOutput("orig_plot_list_ui"))
          )
        ),
        div(
          class = "main-content",
          div(
            class = "content-header",
            div(
              div(class = "content-title", textOutput("orig_selected_plot_title")),
              div(
                class = "breadcrumb-line",
                span(class = "crumb-accent", "Original Pipeline"),
                HTML("&nbsp;&nbsp;&gt;&nbsp;&nbsp;"),
                span(textOutput("orig_selected_plot_breadcrumb"))
              )
            )
          ),
          div(class = "plot-panel", plotOutput("orig_plot", height = "760px")),
          div(
            class = "annotation-panel",
            div(class = "annotation-title", "Annotation or Caption"),
            div(class = "annotation-text", textOutput("orig_annotation"))
          )
        )
      )
    } else if (current_page() == "alternative") {
      div(
        class = "app-shell",
        div(
          class = "sidebar",
          div(
            class = "sidebar-card",
            div(class = "side-section-title", "Alternative Pipeline"),
            div(
              class = "plot-select-wrap",
              selectInput(
                "alt_plot_choice",
                NULL,
                choices = names(alternative_plots),
                selected = if (is.null(input$alt_plot_choice)) names(alternative_plots)[1] else input$alt_plot_choice,
                width = "100%"
              )
            ),
            div(class = "plot-list-wrap", uiOutput("alt_plot_list_ui"))
          )
        ),
        div(
          class = "main-content",
          div(
            class = "content-header",
            div(
              div(class = "content-title", textOutput("alt_selected_plot_title")),
              div(
                class = "breadcrumb-line",
                span(class = "crumb-accent", "Alternative Pipeline"),
                HTML("&nbsp;&nbsp;&gt;&nbsp;&nbsp;"),
                span(textOutput("alt_selected_plot_breadcrumb"))
              )
            )
          ),
          div(class = "plot-panel", plotOutput("alt_plot", height = "760px")),
          div(
            class = "annotation-panel",
            div(class = "annotation-title", "Annotation or Caption"),
            div(class = "annotation-text", textOutput("alt_annotation"))
          )
        )
      )
    } else if (current_page() == "interactive") {
      div(
        class = "interactive-shell",
        div(
          class = "interactive-card",
          fluidRow(
            column(
              3,
              selectInput(
                "interactive_plot_choice",
                "Choose Plot",
                choices = c(
                  "PCA Overlap" = "pca",
                  "Organoid Cluster Comparison" = "org",
                  "Original Figure 3D" = "orig3d",
                  "Alternative Figure 3D" = "alt3d",
                  "Non Coding Figure 3D" = "nc3d"
                )
              ),
              sliderInput("interactive_plot_width", "Plot Width Percent", min = 50, max = 100, value = 100, step = 5),
              sliderInput("interactive_plot_height", "Plot Height Pixels", min = 500, max = 1200, value = 700, step = 50),
              numericInput("n_rows", "Rows of cell assignment table", value = 20, min = 5, max = nrow(cell_data), step = 5),
              checkboxInput("show_table", "Show cell assignment table", TRUE)
            ),
            column(
              9,
              uiOutput("interactive_plot_ui"),
              div(
                class = "axis-card",
                h4("X and Y Axis Table"),
                uiOutput("interactive_axis_table_ui")
              )
            )
          )
        ),
        conditionalPanel(
          condition = "input.show_table == true",
          div(class = "interactive-card", h3("Cell Assignment Table"), tableOutput("cell_table"))
        )
      )
    } else if (current_page() == "gene") {
      div(
        class = "gene-layout",
        div(
          class = "content-header",
          div(
            div(class = "content-title", "Gene Query"),
            div(class = "breadcrumb-line", span(class = "crumb-accent", "Expression Matrix"), HTML("&nbsp;&nbsp;&gt;&nbsp;&nbsp;"), "Individual gene dropdown")
          )
        ),
        div(
          class = "gene-page-grid",
          div(
            class = "gene-control-card",
            h4("Select Gene"),
            p(class = "small-muted", "Choose a gene from the expression matrix dropdown."),
            uiOutput("gene_dropdown_ui"),
            br(),
            p(class = "small-muted", textOutput("gene_matrix_status"))
          ),
          div(
            div(
              class = "gene-result-card",
              h3(textOutput("gene_result_title")),
              tableOutput("gene_summary_table")
            ),
            div(
              class = "gene-result-card",
              h3("Expression Distribution"),
              plotOutput("gene_expression_plot", height = "520px")
            ),
            div(
              class = "gene-result-card",
              h3("Preview of Expression Values"),
              tableOutput("gene_preview_table")
            )
          )
        )
      )
    } else {
      div(
        class = "interactive-shell",
        div(
          class = "conclusion-box",
          h2("Conclusion"),
          div(style = "white-space: pre-wrap; line-height: 1.8; color: #4a5a75;", conclusion_text)
        )
      )
    }
  })
  
  output$orig_selected_plot_title <- renderText({
    req(input$orig_plot_choice)
    input$orig_plot_choice
  })
  
  output$orig_selected_plot_breadcrumb <- renderText({
    req(input$orig_plot_choice)
    input$orig_plot_choice
  })
  
  output$alt_selected_plot_title <- renderText({
    req(input$alt_plot_choice)
    input$alt_plot_choice
  })
  
  output$alt_selected_plot_breadcrumb <- renderText({
    req(input$alt_plot_choice)
    input$alt_plot_choice
  })
  
  output$orig_plot <- renderPlot({
    req(input$orig_plot_choice)
    
    files <- original_plots[[input$orig_plot_choice]]
    env <- plot_cache()
    key <- paste0("original::", input$orig_plot_choice)
    
    if (!exists(key, envir = env, inherits = FALSE)) {
      found <- files[file.exists(files)]
      if (length(found) == 0) stop(paste("Missing plot file for", input$orig_plot_choice))
      assign(key, readRDS(found[1]), envir = env)
      plot_cache(env)
    }
    
    obj <- get(key, envir = env, inherits = FALSE)
    
    if ("recordedplot" %in% class(obj)) {
      replayPlot(obj)
    } else {
      print(obj)
    }
  }, res = 110, width = function() session$clientData$output_orig_plot_width, height = 760)
  
  output$alt_plot <- renderPlot({
    req(input$alt_plot_choice)
    
    files <- alternative_plots[[input$alt_plot_choice]]
    env <- plot_cache()
    key <- paste0("alternative::", input$alt_plot_choice)
    
    if (!exists(key, envir = env, inherits = FALSE)) {
      found <- files[file.exists(files)]
      if (length(found) == 0) stop(paste("Missing plot file for", input$alt_plot_choice))
      assign(key, readRDS(found[1]), envir = env)
      plot_cache(env)
    }
    
    obj <- get(key, envir = env, inherits = FALSE)
    
    if ("recordedplot" %in% class(obj)) {
      replayPlot(obj)
    } else {
      print(obj)
    }
  }, res = 110, width = function() session$clientData$output_alt_plot_width, height = 760)
  
  output$orig_annotation <- renderText({
    req(input$orig_plot_choice)
    if (is.null(plot_notes[[input$orig_plot_choice]])) "" else plot_notes[[input$orig_plot_choice]]
  })
  
  output$alt_annotation <- renderText({
    req(input$alt_plot_choice)
    if (is.null(plot_notes[[input$alt_plot_choice]])) "" else plot_notes[[input$alt_plot_choice]]
  })
  
  output$interactive_plot_ui <- renderUI({
    req(input$interactive_plot_width, input$interactive_plot_height)
    div(
      style = paste0("width:", input$interactive_plot_width, "%; margin: 0 auto;"),
      div(
        class = "plot-panel",
        plotOutput("interactive_plot", height = paste0(input$interactive_plot_height, "px"))
      )
    )
  })
  
  output$interactive_plot <- renderPlot({
    req(input$interactive_plot_choice)
    
    files <- switch(
      input$interactive_plot_choice,
      "pca" = c("PCA_overlap_plot.rds"),
      "org" = c("Organoidclustercomparison.rds"),
      "orig3d" = c("Fig3D_plot.rds"),
      "alt3d" = c("Alt_Fig3D.rds", "Alt_Fig3D - Copy.rds"),
      "nc3d" = c("AltnoncodingFig3D.rds", "AltnoncodingFig3D - Copy.rds")
    )
    
    key <- paste0("interactive::", input$interactive_plot_choice)
    env <- plot_cache()
    
    if (!exists(key, envir = env, inherits = FALSE)) {
      found <- files[file.exists(files)]
      if (length(found) == 0) stop("Missing interactive plot file")
      assign(key, readRDS(found[1]), envir = env)
      plot_cache(env)
    }
    
    obj <- get(key, envir = env, inherits = FALSE)
    
    if ("recordedplot" %in% class(obj)) {
      replayPlot(obj)
    } else {
      print(obj)
    }
  }, res = 110, width = function() session$clientData$output_interactive_plot_width, height = function() input$interactive_plot_height)
  
  output$interactive_axis_table_ui <- renderUI({
    req(input$interactive_plot_choice)
    df <- interactive_axis_info[[input$interactive_plot_choice]]
    
    tags$div(
      class = "axis-table-wrap",
      tags$table(
        class = "axis-table",
        tags$thead(
          tags$tr(
            tags$th("Axis"),
            tags$th("Description")
          )
        ),
        tags$tbody(
          lapply(seq_len(nrow(df)), function(i) {
            tags$tr(
              tags$td(class = "axis-col", df$Axis[i]),
              tags$td(class = "desc-col", df$Description[i])
            )
          })
        )
      )
    )
  })
  
  output$cell_table <- renderTable({
    req(input$n_rows)
    head(cell_data, input$n_rows)
  }, striped = TRUE, bordered = TRUE, spacing = "s")
  
  output$gene_matrix_status <- renderText({
    paste("Matrix file:", gene_matrix_file)
  })
  
  output$gene_dropdown_ui <- renderUI({
    obj <- gene_data_cache()
    
    if (is.null(obj)) {
      if (!file.exists(gene_matrix_file)) {
        return(paste("Gene matrix file not found:", gene_matrix_file))
      }
      
      obj <- readRDS(gene_matrix_file)
      gene_data_cache(obj)
    }
    
    if (is.matrix(obj) || inherits(obj, "Matrix")) {
      if (!is.null(rownames(obj))) {
        genes <- rownames(obj)
      } else {
        genes <- colnames(obj)
      }
    } else if (is.data.frame(obj)) {
      genes <- colnames(obj)
    } else {
      genes <- character(0)
    }
    
    genes <- sort(unique(genes))
    genes <- genes[!is.na(genes)]
    genes <- genes[genes != ""]
    
    selectInput(
      "selected_gene",
      "Select Gene",
      choices = genes,
      selected = if ("SOX2" %in% genes) "SOX2" else genes[1],
      width = "100%"
    )
  })
  
  selected_gene_result <- reactive({
    req(input$selected_gene)
    
    gene_name <- input$selected_gene
    obj <- gene_data_cache()
    
    if (is.null(obj)) {
      req(file.exists(gene_matrix_file))
      obj <- readRDS(gene_matrix_file)
      gene_data_cache(obj)
    }
    
    expression_values <- NULL
    
    if (is.matrix(obj) || inherits(obj, "Matrix")) {
      if (!is.null(rownames(obj)) && gene_name %in% rownames(obj)) {
        expression_values <- as.numeric(obj[gene_name, ])
      } else if (!is.null(colnames(obj)) && gene_name %in% colnames(obj)) {
        expression_values <- as.numeric(obj[, gene_name])
      }
    } else if (is.data.frame(obj)) {
      if (gene_name %in% colnames(obj)) {
        expression_values <- as.numeric(obj[[gene_name]])
      } else if (!is.null(rownames(obj)) && gene_name %in% rownames(obj)) {
        expression_values <- as.numeric(obj[gene_name, ])
      }
    }
    
    if (is.null(expression_values)) {
      list(found = FALSE, gene = gene_name)
    } else {
      expression_values <- expression_values[!is.na(expression_values)]
      list(found = TRUE, gene = gene_name, values = expression_values)
    }
  })
  
  output$gene_result_title <- renderText({
    res <- selected_gene_result()
    
    if (!isTRUE(res$found)) {
      paste("Gene not found:", res$gene)
    } else {
      paste("Gene:", res$gene)
    }
  })
  
  output$gene_summary_table <- renderTable({
    res <- selected_gene_result()
    
    if (!isTRUE(res$found)) {
      return(data.frame(Message = paste("The gene", res$gene, "was not found in the expression matrix.")))
    }
    
    values <- res$values
    
    data.frame(
      Metric = c("Gene", "Number of cells or samples", "Detected values above zero", "Mean expression", "Median expression", "Maximum expression"),
      Value = c(
        res$gene,
        length(values),
        sum(values > 0),
        round(mean(values), 4),
        round(median(values), 4),
        round(max(values), 4)
      ),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")
  
  output$gene_expression_plot <- renderPlot({
    res <- selected_gene_result()
    
    if (!isTRUE(res$found)) {
      plot.new()
      text(0.5, 0.5, "No valid gene selected.")
      return()
    }
    
    values <- res$values
    
    par(mar = c(5, 5, 4, 2))
    hist(
      values,
      breaks = 40,
      main = paste("Expression Distribution for", res$gene),
      xlab = "Expression value",
      ylab = "Number of cells or samples"
    )
  }, res = 110)
  
  output$gene_preview_table <- renderTable({
    res <- selected_gene_result()
    
    if (!isTRUE(res$found)) {
      return(data.frame(Message = "No expression values to show."))
    }
    
    values <- res$values
    n_preview <- min(20, length(values))
    
    data.frame(
      Index = seq_len(n_preview),
      Expression = round(values[seq_len(n_preview)], 4)
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")
  
  observe({
    session$sendCustomMessage("toggleNavActive", list(page = current_page()))
  })
}

js <- HTML("
Shiny.addCustomMessageHandler('toggleNavActive', function(message) {
  var original = document.getElementById('nav_original');
  var alternative = document.getElementById('nav_alternative');
  var interactive = document.getElementById('nav_interactive');
  var gene = document.getElementById('nav_gene');
  var conclusion = document.getElementById('nav_conclusion');

  if (original) original.classList.remove('active-tab');
  if (alternative) alternative.classList.remove('active-tab');
  if (interactive) interactive.classList.remove('active-tab');
  if (gene) gene.classList.remove('active-tab');
  if (conclusion) conclusion.classList.remove('active-tab');

  if (message.page === 'original' && original) original.classList.add('active-tab');
  if (message.page === 'alternative' && alternative) alternative.classList.add('active-tab');
  if (message.page === 'interactive' && interactive) interactive.classList.add('active-tab');
  if (message.page === 'gene' && gene) gene.classList.add('active-tab');
  if (message.page === 'conclusion' && conclusion) conclusion.classList.add('active-tab');
});
")

ui <- tagList(ui_base, tags$script(js))

shinyApp(ui, server)
