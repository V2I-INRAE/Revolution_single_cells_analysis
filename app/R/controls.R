# Controls are a Shiny module; the data and plotting modules do not depend on it.
explorer_group_choices <- function(data, column) {
  groups <- levels(data$cells[[column]])
  counts <- table(data$cells[[column]])
  labels <- groups
  list(values = groups, names = lapply(seq_along(groups), function(i) {
    shiny::tags$span(class = "group-choice",
      shiny::tags$span(class = "group-swatch",
        style = paste0("background:", data$colors[[column]][[groups[i]]])),
      shiny::tags$span(class = "group-name", gsub("_", " ", labels[i])),
      shiny::tags$span(class = "group-count",
        format(counts[[groups[i]]], big.mark = ",")))
  }))
}

explorer_controls_ui <- function(id, views) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::selectInput(ns("view"), "Dataset / view",
      choices = stats::setNames(names(views),
        vapply(views, function(view) view$label, character(1)))),
    shiny::radioButtons(ns("mode"), "Color by",
      choices = c("Clusters" = "cluster", "Gene expression" = "gene"),
      selected = "cluster", inline = TRUE),
    shiny::conditionalPanel("input.mode === 'gene'", ns = ns,
      shiny::selectizeInput(ns("gene"), "Gene", choices = NULL,
        options = list(placeholder = "Search a symbol or Ensembl ID",
          maxOptions = 40L)),
      shiny::checkboxInput(ns("restrict"), "Color checked groups only", FALSE),
      shiny::tags$p(class = "control-help",
        "Expression uses existing RNA LogNormalize values. Zero is low purple.")),
    shiny::conditionalPanel("input.view === 'parent'", ns = ns,
      shiny::selectInput(ns("grouping"), "Group cells by",
        choices = c("Cluster number" = "cluster"))),
    shiny::div(class = "selection-actions",
      shiny::actionButton(ns("all"), "Select all", class = "btn-sm"),
      shiny::actionButton(ns("clear"), "Clear", class = "btn-sm")),
    shiny::checkboxGroupInput(ns("groups"), "Highlight groups",
      choices = character(0)),
    shiny::tags$p(class = "control-help",
      "Unchecked groups remain black. Changing the view or grouping clears the selection.")
  )
}

explorer_controls_server <- function(id, data) {
  shiny::moduleServer(id, function(input, output, session) {
    grouping <- shiny::reactive({
      if (identical(input$view, "parent")) input$grouping else "cluster"
    })
    shiny::observeEvent(data(), {
      current <- data()
      selected <- shiny::isolate(input$gene)
      if (!shiny::isTruthy(selected) || !selected %in% current$genes$id) selected <- ""
      shiny::updateSelectizeInput(session, "gene",
        choices = stats::setNames(current$genes$id, current$genes$label),
        selected = selected, server = TRUE)
    })
    shiny::observeEvent(list(data(), grouping()), {
      choices <- explorer_group_choices(data(), grouping())
      shiny::freezeReactiveValue(input, "groups")
      shiny::updateCheckboxGroupInput(session, "groups",
        choiceNames = choices$names, choiceValues = choices$values,
        selected = character(0))
    })
    shiny::observeEvent(input$all, {
      shiny::updateCheckboxGroupInput(session, "groups",
        selected = levels(data()$cells[[grouping()]]))
    })
    shiny::observeEvent(input$clear, {
      shiny::updateCheckboxGroupInput(session, "groups", selected = character(0))
    })
    list(view = shiny::reactive(input$view),
      mode = shiny::reactive(input$mode), gene = shiny::reactive(input$gene),
      grouping = grouping,
      groups = shiny::reactive(if (is.null(input$groups)) character(0) else input$groups),
      restrict = shiny::reactive(isTRUE(input$restrict)))
  })
}
