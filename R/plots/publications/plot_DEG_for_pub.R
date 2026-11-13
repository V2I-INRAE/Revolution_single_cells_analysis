# Combine saved per-cluster DEG tables of one contrast into one figure; no refit.
# Called by R/plots/publications/main.R. Final-size PDF: 170 mm wide for up to five clusters; place at 100%.
suppressPackageStartupMessages({
  library(ggplot2)
  library(plotthis)
})
source("R/utils/save_plots.R")
source("R/plots/publications/publication_style.R")

plot_DEG_for_pub <- function(
  input_dir,
  output_dir,
  contrast,
  population,
  resolution
) {
  fdr <- deg_config$deg_fdr
  lfc <- deg_config$deg_logfc
  seed <- deg_config$seed
  folders <- sort(Sys.glob(file.path(input_dir, "cluster_[0-9]*")))
  ids <- as.integer(sub("^cluster_", "", basename(folders)))
  paths <- file.path(folders, contrast, paste0(contrast, "_all.csv"))
  tested <- file.exists(paths)
  # Unsupported comparisons must not be displayed as tested results.
  if (!any(tested)) {
    message("No cluster tested ", contrast, "; no figure.")
    return(invisible(NULL))
  }
  paths <- paths[tested]
  ids <- ids[tested]
  dat <- do.call(
    rbind,
    lapply(seq_along(paths), function(i) {
      d <- read.csv(paths[i], stringsAsFactors = FALSE)
      d$cluster <- ids[i]
      d
    })
  )
  dat$class <- ifelse(
    dat$FDR > fdr,
    "Not significant",
    ifelse(
      abs(dat$logFC) < lfc,
      "Significant, small effect",
      ifelse(dat$logFC < 0, "Down", "Up")
    )
  )
  class_colours <- c(
    "Below threshold" = "#CCCCCC",
    "Down" = "#88699B",
    "Up" = "#44AA99"
  )
  dat$plot_class <- factor(
    ifelse(dat$class %in% c("Down", "Up"), dat$class, "Below threshold"),
    levels = names(class_colours)
  )
  # Ensembl IDs and LOC placeholders are plotted but not labelled.
  named <- nzchar(dat$display_gene) &
    !grepl("^(ENS[A-Z]*G[0-9]+([.][0-9]+)?|LOC[0-9]+)$", dat$display_gene)
  dat$annotate <- FALSE
  for (cl in ids) {
    for (direction in c("Down", "Up")) {
      pool <- which(dat$cluster == cl & dat$class == direction & named)
      pool <- pool[order(
        dat$FDR[pool],
        dat$PValue[pool],
        -abs(dat$logFC[pool]),
        dat$gene[pool]
      )]
      dat$annotate[head(pool, 5)] <- TRUE
    }
  }
  set.seed(seed)
  dat$x <- match(dat$cluster, ids) - 0.32 + runif(nrow(dat), -0.10, 0.10)
  dat <- dat[order(dat$plot_class), ] # Genes passing both thresholds drawn on top.
  strip <- data.frame(
    x = seq_along(ids),
    label = paste("Cluster", ids),
    colour = viz_palette[ids]
  )
  strip$text <- ifelse(
    colSums(grDevices::col2rgb(strip$colour) * c(0.299, 0.587, 0.114)) > 160,
    "black",
    "white"
  )
  p <- ScatterPlot(
    dat,
    x = "x",
    y = "logFC",
    color_by = "plot_class",
    palcolor = class_colours,
    shape = 16,
    size_by = 0.65,
    alpha = 1,
    border_size = 0,
    seed = seed,
    theme = "theme_classic",
    legend.position = "bottom",
    aspect.ratio = NULL,
    combine = FALSE
  )[[1]] +
    geom_hline(
      yintercept = c(-lfc, lfc),
      linetype = "dashed",
      colour = "grey60",
      linewidth = 0.25
    ) +
    geom_rect(
      data = strip,
      aes(
        xmin = x - 0.5,
        xmax = x + 0.5,
        ymin = -0.14,
        ymax = 0.14,
        fill = I(colour)
      ),
      inherit.aes = FALSE
    ) +
    geom_text(
      data = strip,
      aes(x, 0, label = label, colour = I(text)),
      inherit.aes = FALSE,
      show.legend = FALSE,
      size = 11 / ggplot2::.pt
    ) +
    scale_x_continuous(limits = c(0.45, length(ids) + 0.55), breaks = NULL) +
    scale_y_continuous(
      limits = range(dat$logFC) + c(-1.5, 1.5),
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0(population, ": ", gsub("_vs_", " vs ", contrast)),
      x = NULL,
      subtitle = paste0(
        "Saved DEG: resolution ",
        resolution,
        "; clusters ",
        paste(ids, collapse = ", ")
      ),
      y = paste0("Log2 fold change (", gsub("_vs_", " / ", contrast), ")"),
      colour = NULL,
      caption = deg_wrap(
        sprintf(
          "Dashed guides: log2FC +/-%g. FDR <= %g within each cluster; no global FDR adjustment.",
          lfc,
          fdr
        ),
        width = 80
      )
    ) +
    deg_theme() +
    theme(
      axis.line.x = element_blank(),
      legend.position = "bottom",
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 12)
    ) +
    guides(
      colour = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(size = 2)
      )
    )
  for (i in seq_along(ids)) {
    for (direction in c("Down", "Up")) {
      lab <- dat[
        dat$cluster == ids[i] & dat$annotate & dat$class == direction,
      ]
      p <- p +
        ggrepel::geom_text_repel(
          data = lab,
          aes(x, logFC, label = display_gene),
          inherit.aes = FALSE,
          show.legend = FALSE,
          seed = seed,
          size = 9.5 / ggplot2::.pt,
          max.overlaps = Inf,
          xlim = c(i - 0.5, i + 0.75),
          direction = "y",
          hjust = 0,
          nudge_x = i - 0.18 - lab$x,
          nudge_y = sign(lab$logFC) * 0.6,
          ylim = if (direction == "Up") {
            c(0.55, max(dat$logFC) + 1)
          } else {
            c(min(dat$logFC) - 1, -0.55)
          },
          box.padding = 0.35,
          force = 10,
          force_pull = 0,
          segment.colour = "grey40",
          segment.size = 0.2,
          min.segment.length = 0,
          max.time = Inf,
          max.iter = 100000
        )
    }
  }
  p <- ggrastr::rasterise(p, layers = "Point", dpi = 600)
  out <- file.path(output_dir, paste0("DEG_", contrast, ".pdf"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    out,
    p,
    device = grDevices::cairo_pdf,
    width = max(170, length(ids) * 34),
    height = 170,
    units = "mm",
    bg = "white"
  )
  message("Saved: ", out)
  invisible(out)
}
