# One GSEA collection and contrast across all tested clusters; no refit.
# Reads <deg_dir>/cluster_XX/<contrast>/gsea/<collection>/<contrast>_gsea.csv.
# Labels: whitelisted significant terms first (exact names in the CSV; edit it), then the
# lowest-FDR significant terms (italic). Reactome/GO:BP labels use main pathways only.
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})
source("R/plots/publications/publication_style.R")

gsea_overview_labels <- c(
  hallmark = "Hallmark",
  reactome = "Reactome",
  gobp = "GO:BP",
  kegg = "KEGG"
)

plot_gsea_overview <- function(
  deg_dir,
  output_dir,
  collection,
  contrast,
  population,
  n_labels = 8,
  whitelist = "R/plots/publications/gsea_label_whitelist.csv"
) {
  fdr <- deg_config$pathway_fdr
  seed <- deg_config$seed
  groups <- strsplit(contrast, "_vs_", fixed = TRUE)[[1]]
  collection_label <- gsea_overview_labels[[collection]]
  folders <- sort(Sys.glob(file.path(deg_dir, "cluster_[0-9]*")))
  paths <- file.path(
    folders,
    contrast,
    "gsea",
    collection,
    paste0(contrast, "_gsea.csv")
  )
  # Skipped or blocked cluster/contrast pairs have no table and are not shown.
  folders <- folders[file.exists(paths)]
  paths <- paths[file.exists(paths)]
  if (!length(paths)) {
    return(invisible(NULL))
  }
  ids <- as.integer(sub("^cluster_", "", basename(folders)))
  dat <- do.call(
    rbind,
    lapply(seq_along(paths), function(i) {
      x <- read.csv(paths[i])
      x$cluster <- ids[i]
      x[c("pathway", "padj", "NES", "size", "main_pathway", "cluster")]
    })
  )
  dat$significant <- dat$padj < fdr
  allowed <- read.csv(whitelist)
  dat$whitelisted <- dat$pathway %in%
    allowed$pathway[allowed$collection == collection]
  # Cleaned term name, as in the dot plots; GO descriptions are definitions, not names.
  dat$label <- gsub(
    "_",
    " ",
    sub("^(HALLMARK|REACTOME|GOBP|KEGG)_", "", dat$pathway)
  )
  dat$label <- ifelse(
    nchar(dat$label) > 50,
    paste0(substr(dat$label, 1, 47), "..."),
    dat$label
  )
  dat$label <- vapply(
    dat$label,
    function(x) paste(strwrap(x, width = 22), collapse = "\n"),
    character(1)
  )
  # Per cluster and direction: whitelisted first, then lowest FDR; main pathways only
  # (fgsea::collapsePathways groups the other significant terms under them).
  dat$annotate <- FALSE
  dat$face <- "plain"
  for (cl in ids) {
    for (positive in c(TRUE, FALSE)) {
      pool <- which(
        dat$cluster == cl &
          dat$significant &
          dat$main_pathway %in% TRUE &
          (dat$NES > 0) == positive
      )
      pool <- pool[order(
        !dat$whitelisted[pool],
        dat$padj[pool],
        -abs(dat$NES[pool]),
        dat$pathway[pool]
      )]
      chosen <- head(pool, n_labels)
      dat$annotate[chosen] <- TRUE
      dat$face[chosen] <- ifelse(dat$whitelisted[chosen], "plain", "italic")
    }
  }
  set.seed(seed)
  dat$x <- match(dat$cluster, ids) + runif(nrow(dat), -0.28, 0.28)
  dat <- dat[order(dat$significant), ] # Significant points drawn on top.
  colour <- viz_palette[[match(collection, deg_config$collections)]]
  # Room above/below the points for the label stacks, only on sides that have labels.
  span <- max(abs(dat$NES))
  room <- function(up) {
    n <- max(0, table(dat$cluster[dat$annotate & (dat$NES > 0) == up]))
    if (n == 0) 0.05 * span else span * (0.1 + 0.16 * n)
  }
  upper <- max(dat$NES, 0) + room(TRUE)
  lower <- min(dat$NES, 0) - room(FALSE)
  # Labels are stacked above the highest (or below the lowest) point of their cluster.
  label_layers <- unlist(
    lapply(seq_along(ids), function(i) {
      lapply(c(TRUE, FALSE), function(positive) {
        lab <- dat[
          dat$annotate & dat$cluster == ids[i] & (dat$NES > 0) == positive,
        ]
        in_cluster <- dat$NES[dat$cluster == ids[i]]
        target <- if (positive) {
          max(in_cluster) + 0.04 * span
        } else {
          min(in_cluster) - 0.04 * span
        }
        ggrepel::geom_text_repel(
          data = lab,
          aes(label = label, fontface = face),
          size = 6.5 / ggplot2::.pt,
          lineheight = 0.85,
          seed = seed,
          max.overlaps = Inf,
          box.padding = 0.35,
          min.segment.length = 0,
          segment.size = 0.2,
          segment.colour = "grey40",
          force = 10,
          max.time = 10,
          max.iter = 100000,
          direction = "both",
          xlim = c(i - 0.45, i + 0.45),
          ylim = if (positive) c(target, upper) else c(lower, target),
          nudge_x = (i - lab$x) * 0.7,
          nudge_y = target - lab$NES
        )
      })
    }),
    recursive = FALSE
  )
  tested <- table(factor(dat$cluster, levels = ids))
  p <- ggplot(dat, aes(x, NES)) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.3
    ) +
    geom_vline(
      xintercept = seq_along(ids)[-length(ids)] + 0.5,
      colour = "grey70",
      linewidth = 0.25
    ) +
    geom_point(
      data = dat[!dat$significant, ],
      aes(size = size),
      shape = 21,
      fill = colour,
      colour = "transparent",
      alpha = 0.25
    ) +
    geom_point(
      data = dat[dat$significant, ],
      aes(size = size),
      shape = 21,
      fill = colour,
      colour = "black",
      stroke = 0.3,
      alpha = 0.9
    ) +
    label_layers +
    scale_size_continuous(range = c(0.8, 4), name = "Matched genes") +
    scale_x_continuous(
      breaks = seq_along(ids),
      labels = paste("Cluster", ids),
      limits = c(0.5, length(ids) + 0.5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(lower, upper), expand = c(0, 0)) +
    labs(
      title = paste0(
        population,
        ": ",
        collection_label,
        " GSEA, ",
        gsub("_vs_", " vs ", contrast)
      ),
      subtitle = paste0(
        groups[1],
        " minus ",
        groups[2],
        " | solid: pathway FDR < ",
        fdr,
        " | up to ",
        n_labels,
        " labels per cluster and direction"
      ),
      x = NULL,
      y = paste0(
        "NES\n",
        groups[2],
        " higher (down) | ",
        groups[1],
        " higher (up)"
      ),
      caption = deg_wrap(
        paste0(
          "Plain labels were pre-selected by name (inflammation, stress, cell death); italic labels are the lowest-FDR significant terms. ",
          "Only main pathways are labelled: fgsea::collapsePathways groups significant terms explained by a more significant pathway under it. ",
          "BH FDR within each cluster, collection and contrast; no adjustment across clusters.\n",
          "Pathways tested: ",
          paste(paste0("cluster ", ids, " ", tested), collapse = ", "),
          "."
        ),
        width = 100
      )
    ) +
    deg_theme() +
    theme(
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(size = 10, face = "bold")
    )
  save_deg_figure(
    p,
    output_dir,
    paste0("GSEA_", collection, "_", contrast),
    width = max(170, length(ids) * 34),
    height = 190
  )
  invisible(p)
}
