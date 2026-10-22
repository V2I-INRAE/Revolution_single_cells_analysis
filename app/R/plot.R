# Plotly traces use the saved embedding, with no dimensionality reduction here.
explorer_ranges <- function(data) {
  lapply(data$cells[c("umap_1", "umap_2")], function(x) {
    bounds <- range(x)
    bounds + c(-1, 1) * diff(bounds) * 0.04
  })
}

explorer_foreground <- function(data, grouping, groups, gene = NULL,
  restrict = FALSE) {
  selected <- data$cells[[grouping]] %in% groups
  indices <- which(if (is.null(gene) || restrict) selected else
    rep(TRUE, nrow(data$cells)))
  marker <- list(size = 2, opacity = 1, showscale = FALSE, color = "#000000")
  hover <- "%{text}<extra></extra>"
  if (is.null(gene)) {
    marker$color <- unname(data$colors[[grouping]][
      as.character(data$cells[[grouping]][indices])])
    text <- paste("Group", data$cells[[grouping]][indices])
  } else {
    # High expression is drawn last, as in the existing ordered FeaturePlots.
    indices <- indices[order(gene$values[indices], method = "radix")]
    palette <- viridisLite::viridis(256)
    marker$color <- gene$values[indices]
    marker$colorscale <- lapply(seq_along(palette), function(i) {
      list((i - 1L) / (length(palette) - 1L), palette[i])
    })
    marker$cmin <- 0
    marker$cmax <- if (gene$maximum > 0) gene$maximum else 1
    marker$showscale <- TRUE
    marker$colorbar <- list(title = list(text = "RNA<br>LogNormalize"),
      thickness = 14, len = 0.65, outlinewidth = 0)
    text <- rep(gene$label, length(indices))
    hover <- "%{text}<br>Normalized RNA: %{marker.color:.3f}<extra></extra>"
  }
  list(x = data$cells$umap_1[indices], y = data$cells$umap_2[indices],
    marker = marker, text = text, hovertemplate = hover)
}

explorer_plot <- function(data) {
  ranges <- explorer_ranges(data)
  plotly::plot_ly(source = "umap") |>
    plotly::add_trace(type = "scattergl", mode = "markers",
      x = data$cells$umap_1, y = data$cells$umap_2,
      marker = list(color = "#000000", size = 2, opacity = 1),
      hoverinfo = "skip", showlegend = FALSE, name = "All cells") |>
    # Keep a second trace even before any cells are highlighted.
    plotly::add_trace(type = "scattergl", mode = "markers",
      x = NA_real_, y = NA_real_,
      marker = list(color = "#000000", size = 2, opacity = 1),
      hoverinfo = "text", showlegend = FALSE, name = "Highlighted cells") |>
    # Different embeddings reset the axes; colors within a view preserve zoom.
    plotly::layout(dragmode = "pan", uirevision = data$id,
      paper_bgcolor = "white", plot_bgcolor = "white",
      # Reserve colorbar space so toggling expression cannot change axis ranges.
      margin = list(l = 55, r = 140, t = 15, b = 50, autoexpand = FALSE),
      font = list(family = "Arial, sans-serif", color = "#334155"),
      xaxis = list(title = "UMAP 1", range = ranges$umap_1,
        zeroline = FALSE, showgrid = FALSE),
      yaxis = list(title = "UMAP 2", range = ranges$umap_2,
        zeroline = FALSE, showgrid = FALSE, scaleanchor = "x", scaleratio = 1)) |>
    plotly::config(scrollZoom = TRUE, displaylogo = FALSE,
      modeBarButtonsToRemove = c("select2d", "lasso2d", "autoScale2d")) |>
    htmlwidgets::onRender("function(el, x, view) {
      // Remove only this view's handler, including a pending animation frame.
      if (el._explorerZoomCleanup) el._explorerZoomCleanup();
      el._explorerZoomCleanup = null;
      if (!view.parent) {
        let frame = null;
        const resize = function() {
          if (frame !== null) return;
          frame = requestAnimationFrame(function() {
            frame = null;
            const axes = el._fullLayout;
            const zoom = Math.max(1,
              (view.ranges.umap_1[1] - view.ranges.umap_1[0]) /
                (axes.xaxis.range[1] - axes.xaxis.range[0]),
              (view.ranges.umap_2[1] - view.ranges.umap_2[0]) /
                (axes.yaxis.range[1] - axes.yaxis.range[0]));
            const size = Math.min(8, 2 * Math.sqrt(zoom));
            // Gene/highlight updates replace the foreground marker; restore its
            // zoom size without another redraw when both sizes already match.
            if (el.data.some(trace => trace.marker.size !== size)) {
              Plotly.restyle(el, {'marker.size': size}, [0, 1]);
            }
          });
        };
        el.on('plotly_afterplot', resize);
        el._explorerZoomCleanup = function() {
          el.removeListener('plotly_afterplot', resize);
          if (frame !== null) cancelAnimationFrame(frame);
        };
      }
      Shiny.setInputValue('umap_ready', x.layout.uirevision, {priority: 'event'});
    }", data = list(parent = data$parent, ranges = ranges))
}

update_explorer_plot <- function(session, foreground) {
  plotly::plotlyProxy("umap", session) |>
    plotly::plotlyProxyInvoke("restyle", list(
      x = list(foreground$x), y = list(foreground$y),
      marker = list(foreground$marker), text = list(foreground$text),
      hovertemplate = foreground$hovertemplate, hoverinfo = "all"), list(1L))
}
