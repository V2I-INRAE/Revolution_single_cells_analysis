# Shiny runs this entry point with app/ as its working directory.
source("config.R", local = TRUE)
source("R/data.R", local = TRUE)
source("R/controls.R", local = TRUE)
source("R/plot.R", local = TRUE)
source("R/export.R", local = TRUE)
views <- explorer_config(normalizePath(".."))[explorer_active_views]
format_count <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

ui <- bslib::page_sidebar(
  title = shiny::div(class = "app-title", "REVOLUTION",
    shiny::span("Cell explorer", class = "app-subtitle")),
  theme = bslib::bs_theme(version = 5, bg = "#F5F6F8", fg = "#18212F",
    primary = "#3A5BA0", base_font = "system-ui"),
  sidebar = bslib::sidebar(title = "Explore cells", width = 345,
    explorer_controls_ui("controls", views)),
  shiny::tags$head(shiny::tags$link(rel = "stylesheet", href = "explorer.css")),
  shiny::uiOutput("dataset", class = "dataset-line"),
  bslib::card(class = "umap-card", full_screen = TRUE,
    bslib::card_header(shiny::div(class = "plot-header",
      shiny::div(shiny::textOutput("heading"),
        shiny::textOutput("summary", container = shiny::tags$small)),
      shiny::div(class = "d-flex gap-2",
        shiny::downloadButton("pdf", "Save PDF", class = "btn-sm",
          title = "600 dpi cell image; editable text, axes and legends."),
        shiny::actionButton("reset", "Reset view", class = "btn-sm")))),
    plotly::plotlyOutput("umap", height = "100%"),
    bslib::card_footer(shiny::uiOutput("notice"),
      shiny::textOutput("source", container = shiny::tags$small)))
)

server <- function(input, output, session) {
  # Only the selected view is opened; sessions do not share mutable selections.
  data <- shiny::reactive({
    id <- controls$view()
    shiny::req(id %in% names(views))
    current <- open_explorer_data(views[[id]])
    current$id <- id
    current
  })
  controls <- explorer_controls_server("controls", data)
  # A view change reads the gene from that view, never from the parent cache.
  gene <- shiny::reactive({
    shiny::req(controls$gene())
    read_explorer_gene(data(), controls$gene())
  })
  current_gene <- shiny::reactive({
    if (identical(controls$mode(), "gene") && shiny::isTruthy(controls$gene()) &&
      controls$gene() %in% data()$genes$id)
      gene() else NULL
  })
  foreground <- shiny::reactive({
    selected <- current_gene()
    groups <- if (controls$mode() == "gene" && is.null(selected)) character(0) else
      controls$groups()
    explorer_foreground(data(), controls$grouping(), groups,
      selected, controls$restrict())
  })
  output$umap <- plotly::renderPlotly(explorer_plot(data()))
  shiny::observe({
    # Wait for this view's black trace before applying its foreground values.
    shiny::req(identical(input$umap_ready, controls$view()))
    update_explorer_plot(session, foreground())
  })
  output$pdf <- shiny::downloadHandler(
    filename = function() {
      current <- data()
      selected <- current_gene()
      feature <- if (is.null(selected)) "clusters" else
        gsub("[^A-Za-z0-9_.-]", "_", selected$id)
      paste0(current$id, "_", current$source$cluster_column, "_", feature, ".pdf")
    },
    content = function(path) {
      current <- data()
      view <- input$umap_view
      shiny::req(identical(view$id, current$id))
      save_explorer_pdf(path, current, foreground(), view, current_gene(),
        current$colors[[controls$grouping()]][controls$groups()])
    },
    contentType = "application/pdf"
  )
  shiny::observeEvent(input$reset, {
    ranges <- explorer_ranges(data())
    plotly::plotlyProxy("umap", session) |>
      plotly::plotlyProxyInvoke("relayout",
        list("xaxis.range" = ranges$umap_1, "yaxis.range" = ranges$umap_2))
  })
  output$dataset <- shiny::renderUI({
    current <- data()
    shiny::tagList(shiny::span(paste(format_count(nrow(current$cells)), "cells")),
      shiny::span(paste(format_count(nrow(current$genes)), "genes")),
      shiny::span(current$label), shiny::span("Saved Harmony UMAP"))
  })
  output$source <- shiny::renderText({
    source <- data()$source
    embedding <- if (!identical(source$embedding_input, source$input))
      paste("| UMAP source:", basename(dirname(source$embedding_input)),
        basename(source$embedding_input), "(barcode matched)") else ""
    paste("Source:", basename(source$input), "|", source$producer, "|",
      source$reduction, "|", source$cluster_column, embedding,
      "| RNA data | No filtering, normalization, or UMAP refit.")
  })
  output$heading <- shiny::renderText({
    selected <- current_gene()
    if (!is.null(selected)) selected$label else if (controls$mode() == "gene")
      "Choose a gene to view its expression" else "Cluster highlighting"
  })
  output$summary <- shiny::renderText({
    current <- data()
    selected <- current_gene()
    checked <- current$cells[[controls$grouping()]] %in% controls$groups()
    if (is.null(selected)) {
      paste(format_count(nrow(current$cells)), "cells shown |",
        format_count(if (controls$mode() == "gene") 0 else sum(checked)), "highlighted")
    } else {
      visible <- if (controls$restrict()) checked else rep(TRUE, nrow(current$cells))
      paste(format_count(sum(visible)), "cells colored |",
        format_count(sum(selected$values[visible] > 0)), "with expression > 0 |",
        "Full-view range: 0 to", format(signif(selected$maximum, 4), trim = TRUE))
    }
  })
  output$notice <- shiny::renderUI({
    selected <- current_gene()
    if (!is.null(selected) && any(selected$unavailable)) {
      shiny::div(class = "availability-notice",
        paste("This feature was absent from the normalized layer for",
          format_count(sum(selected$unavailable)), "cells. Joining layers represents",
          "these values as zero, as in the existing FeaturePlots."))
    }
  })
}

shiny::shinyApp(ui, server)
