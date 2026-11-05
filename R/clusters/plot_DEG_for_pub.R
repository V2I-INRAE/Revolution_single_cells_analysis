# Combine saved per-cluster volcano results; no differential-expression refit.
# Usage: Rscript R/clusters/plot_DEG_for_pub.R <DEG_T10H_directory> [output_dir]
# Or source and call plot_DEG_for_pub(input_dir, output_dir, contrast).
# Final-size PDF: 170 mm wide for up to five clusters; place at 100%, not 20%.
suppressPackageStartupMessages({ library(ggplot2); library(plotthis) })
source("R/utils/save_plots.R")
source("R/DEG/publication_style.R")

plot_DEG_for_pub <- function(input_dir, output_dir = file.path(input_dir, "for_pub"),
  contrast = "P4_vs_N4") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  info_path <- file.path(output_dir, ".INFO")
  info <- c("# Combined-cluster differential-expression plot", "Status: running",
    paste("Started:", Sys.time()), paste("Input (read-only):", normalizePath(input_dir)),
    paste("Producer:", Sys.getenv("SLURM_JOB_ID", unset = "manual/local")),
    paste("Contrast:", contrast, "; source record:", file.path(input_dir, ".INFO")),
    "Method: display saved per-cluster edgeR results; no refitting or pooling of counts.",
    "Thresholds: BH FDR <= 0.05 AND abs(log2FC) >= 0.5; no raw-p cutoff.",
    "FDR families remain per cluster/contrast, not global across the figure.",
    "All tested gene-cluster pairs retained; saved volcano display names/classes reused.",
    "Grey: fails either FDR or fold-change threshold (including significant small effects).",
    "Labels: up to five named genes per cluster/direction; rank FDR, PValue, descending abs(logFC), gene ID.",
    "Ensembl gene IDs and LOC placeholders excluded from labels only; their points remain plotted.",
    "Seed: 1234; gene cloud at cluster-0.32 with jitter +/-0.10; points rasterized at 600 dpi; PDF only.",
    "Positive logFC favours the first condition; negative favours the second.",
    "Typography: plain editable labels (no halo), thin grey connector lines; gene labels 9.5 pt, axes/legend 10 pt, axis titles 11 pt, title 12 pt.",
    "Label layout: plain left-aligned column at cluster-0.18, right of gene cloud; separate signs; vertical repel/nudge +/-0.6, padding0.35, force10, force_pull0; thin grey straight connectors.")
  writeLines(info, info_path)
  done <- FALSE
  on.exit({
    info[info == "Status: running"] <- paste("Status:", if (done) "completed" else "failed/incomplete")
    writeLines(c(info, paste("Finished:", Sys.time())), info_path)
  }, add = TRUE)
  record <- readLines(file.path(input_dir, ".INFO"))
  population_line <- grep("^Population:", record, value = TRUE)
  population <- trimws(sub("^Population: ([^;]+).*", "\\1", population_line))
  cluster_col <- sub(".*cluster column: ([^ ;]+).*", "\\1", population_line)
  resolution <- sub(".*_res[.]", "", cluster_col)

  folders <- list.dirs(input_dir, recursive = FALSE)
  folders <- folders[grepl("^cluster_[0-9]+$", basename(folders))]
  ids <- as.integer(sub("cluster_", "", basename(folders)))
  folders <- folders[order(ids)]; ids <- sort(ids)
  paths <- file.path(folders, "plots", "volcano_data.csv")
  stopifnot(length(ids) > 0, max(ids) <= length(viz_palette), all(file.exists(paths)))
  dat <- do.call(rbind, lapply(seq_along(paths), function(i) {
    d <- read.csv(paths[i], stringsAsFactors = FALSE)
    d <- d[d$contrast == contrast, ]
    stopifnot(nrow(d) > 0, !anyDuplicated(d$gene), !anyNA(d),
      all(is.finite(d$logFC)), all(is.finite(d$FDR)))
    d$cluster <- ids[i]
    d
  }))
  class_colours <- c("Below threshold" = "#CCCCCC", "Down" = "#88699B", "Up" = "#44AA99")
  expected <- ifelse(dat$FDR > 0.05, "Not significant", ifelse(abs(dat$logFC) < 0.5,
    "Significant, small effect", ifelse(dat$logFC < 0, "Down", "Up")))
  stopifnot(identical(dat$class, expected))
  dat$plot_class <- factor(ifelse(dat$class %in% c("Down", "Up"), dat$class, "Below threshold"),
    levels = names(class_colours))
  named <- nzchar(dat$display_gene) & !grepl("^(ENS[A-Z]*G[0-9]+([.][0-9]+)?|LOC[0-9]+)$", dat$display_gene)
  dat$annotate <- FALSE
  for (cl in ids) for (direction in c("Down", "Up")) {
    pool <- which(dat$cluster == cl & dat$class == direction & named)
    pool <- pool[order(dat$FDR[pool], dat$PValue[pool], -abs(dat$logFC[pool]), dat$gene[pool])]
    dat$annotate[head(pool, 5)] <- TRUE
  }
  set.seed(1234)
  dat$x <- match(dat$cluster, ids) - 0.32 + runif(nrow(dat), -0.10, 0.10)
  dat <- dat[order(dat$plot_class), ] # Genes passing both thresholds drawn on top.
  strip <- data.frame(x = seq_along(ids), label = paste("Cluster", ids),
    colour = viz_palette[ids])
  strip$text <- ifelse(colSums(grDevices::col2rgb(strip$colour) * c(0.299, 0.587, 0.114)) > 160,
    "black", "white")
  p <- ScatterPlot(dat, x = "x", y = "logFC", color_by = "plot_class", palcolor = class_colours,
    shape = 16, size_by = 0.65, alpha = 1, border_size = 0, seed = 1234,
    theme = "theme_classic", legend.position = "bottom", aspect.ratio = NULL,
    combine = FALSE)[[1]] +
    geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey60", linewidth = 0.25) +
    geom_rect(data = strip, aes(xmin = x - 0.5, xmax = x + 0.5, ymin = -0.14, ymax = 0.14,
      fill = I(colour)), inherit.aes = FALSE) +
    geom_text(data = strip, aes(x, 0, label = label, colour = I(text)), inherit.aes = FALSE, show.legend = FALSE, size = 11 / ggplot2::.pt) +
    scale_x_continuous(limits = c(0.45, length(ids) + 0.55), breaks = NULL) +
    scale_y_continuous(limits = range(dat$logFC) + c(-1.5, 1.5), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(title = paste0(population, ": ", gsub("_vs_", " vs ", contrast)), x = NULL,
      subtitle = paste0("Saved DEG: resolution ", resolution, "; clusters ", paste(ids, collapse = ", ")),
      y = paste0("Log2 fold change (", gsub("_vs_", " / ", contrast), ")"), colour = NULL,
      caption = deg_wrap("Dashed guides: log2FC +/-0.5. FDR <= 0.05 within each cluster; no global FDR adjustment.", width = 80)) +
    deg_theme() + theme(axis.line.x = element_blank(), legend.position = "bottom",
      axis.text = element_text(size = 10), axis.title = element_text(size = 11),
      legend.text = element_text(size = 10), plot.title = element_text(size = 12)) +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 2)))
  for (i in seq_along(ids)) for (direction in c("Down", "Up")) {
    lab <- dat[dat$cluster == ids[i] & dat$annotate & dat$class == direction, ]
    p <- p + ggrepel::geom_text_repel(data = lab, aes(x, logFC, label = display_gene),
      inherit.aes = FALSE, show.legend = FALSE, seed = 1234, size = 9.5 / ggplot2::.pt, max.overlaps = Inf,
      xlim = c(i - 0.5, i + 0.75), direction = "y", hjust = 0,
      nudge_x = i - 0.18 - lab$x, nudge_y = sign(lab$logFC) * 0.6,
      ylim = if (direction == "Up") c(0.55, max(dat$logFC) + 1) else c(min(dat$logFC) - 1, -0.55),
      box.padding = 0.35, force = 10, force_pull = 0, segment.colour = "grey40",
      segment.size = 0.2, min.segment.length = 0,
      max.time = Inf, max.iter = 100000)
  }
  p <- ggrastr::rasterise(p, layers = "Point", dpi = 600)
  out <- file.path(output_dir, paste0("DEG_all_clusters_", contrast, ".pdf"))
  info <- c(info, paste("Actual inputs:", paste(paths, collapse = "; ")),
    paste("Saved DEG cluster column:", cluster_col, "; clusters plotted:", paste(ids, collapse = ", ")),
    paste("Gene-cluster pairs:", nrow(dat), "; named gene labels:", sum(dat$annotate)),
    paste("Versions: R", getRversion(), "; plotthis", packageVersion("plotthis"),
      "; ggplot2", packageVersion("ggplot2"), "; ggrastr", packageVersion("ggrastr")),
    "Producer analysis: inherited from input .INFO; no changes to the saved model, cohort or adjustment families.",
    paste("PDF dimensions (mm):", max(170, length(ids) * 34), "x 170; place at 100%."),
    capture.output(print(table(dat$cluster, dat$class))), paste("Output:", out))
  writeLines(info, info_path)
  ggsave(out, p, device = grDevices::cairo_pdf, width = max(170, length(ids) * 34),
    height = 170, units = "mm", bg = "white")
  stopifnot(identical(rawToChar(readBin(out, "raw", n = 5)), "%PDF-"))
  info <- c(info, "Validation: saved classes match approved thresholds; complete unique gene rows; nonempty PDF signature checked.",
    "External PDF reopening and visual inspection pending.")
  done <- TRUE
  message("Saved: ", out)
  invisible(list(data = dat, plot = p, pdf = out))
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  if (!length(a) || length(a) > 2L) stop("Usage: Rscript R/clusters/plot_DEG_for_pub.R <DEG_directory> [output_dir]")
  if (length(a) == 1L) plot_DEG_for_pub(a[1]) else plot_DEG_for_pub(a[1], a[2])
}
