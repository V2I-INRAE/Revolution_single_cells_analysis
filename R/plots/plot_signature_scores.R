# Plot saved AUCell or mean log-expression scores; never recompute score values.
# Requested paper-style view: point_mode="cells" draws every scored cell as a
# small black rasterized dot with narrow horizontal jitter; no sample/median dots.
# source("R/plots/plot_signature_scores.R")
# plot_signature_scores("results/signature_scores_parent_new/signature_scores.rds",
#   cluster_col = "harmony_snn_res.0.7", group_by = "cluster",
#   output_root = "results/signature_violins_new")
# For subclusters, supply metadata = readRDS(subset_rds)[[]] and its own cluster_col.
# signatures accepts registry names (both Rooney arms) or exact saved score names.
# Explicit approved five-condition view: group_by="condition",
# condition_grouping="deg_five". REVO31-N8/T8H negative is displayed as N10;
# original sample/time metadata remain unchanged and are exported with plot groups.
source("R/utils/save_plots.R")
source("R/DEG/publication_style.R")
source("R/DEG/signature_conditions.R")

plot_signature_scores <- function(scores_file, cluster_col, output_root, metadata = NULL,
    clusters = NULL, signatures = NULL, group_by = c("cluster", "condition"),
    condition_col = NULL, group_levels = NULL,
    condition_grouping = c("original", "deg_five"),
    point_mode = c("sample_medians", "cells"), point_seed = 1234L,
    display_groups = NULL, comparisons_file = NULL, pooled_mean_bar = FALSE) {
  group_by <- match.arg(group_by)
  point_mode <- match.arg(point_mode)
  if (point_mode == "cells") {
    for (p in c("ggrastr", "withr")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing cell-dot plotting package: ", p)
    if (length(point_seed) != 1L || is.na(point_seed) || point_seed < 0 || point_seed != as.integer(point_seed)) stop("Use a nonnegative integer display-only point_seed.")
  }
  condition_grouping <- match.arg(condition_grouping)
  if (condition_grouping == "deg_five" && (group_by != "condition" || !is.null(condition_col))) {
    stop("deg_five is an explicit condition view; do not combine with condition_col or cluster view.")
  }
  if (file.exists(output_root)) stop("output_root exists; use a new plot directory.")
  if ((!is.null(display_groups) || !is.null(comparisons_file)) &&
      (group_by != "condition" || condition_grouping != "deg_five")) stop("Explicit display selection/comparisons require deg_five condition view.")
  if (pooled_mean_bar && point_mode != "cells") stop("Mean bars require cell-point view without median overlays.")
  x <- readRDS(scores_file)
  score_method <- x$settings$method
  if (is.null(score_method) && isTRUE(x$settings$normAUC) && !is.null(x$settings$auc_max_rank)) score_method <- "AUCell" # existing saved AUCell files
  if (length(score_method) != 1L || !score_method %in% c("AUCell", "mean_log_expression")) stop("Unknown saved scoring method; do not guess expression units.")
  y_label <- if (score_method == "AUCell") "AUCell score (normalized AUC)" else "Mean log-normalized expression"
  scores <- x$scores
  if (anyDuplicated(colnames(scores)) || !identical(colnames(scores), rownames(x$metadata))) stop("Saved score cell IDs invalid.")
  if (is.null(metadata)) metadata <- x$metadata
  required <- c(cluster_col, "sample", "pig", "pressure", "time_point")
  if (!all(required %in% names(metadata)) || anyDuplicated(rownames(metadata)) ||
      anyNA(metadata[, required, drop = FALSE])) stop("Metadata requires unique cell IDs and complete cluster/sample/pig/pressure/time columns.")
  ids <- colnames(scores)[colnames(scores) %in% rownames(metadata)]
  if (!length(ids)) stop("No scored cells match supplied metadata by ID.")
  n_matched <- length(ids)
  md <- metadata[ids, , drop = FALSE]
  # A barcode cannot silently acquire a different specimen identity on a new view.
  for (n in c("sample", "pig", "pressure", "time_point")) {
    if (!identical(as.character(md[[n]]), as.character(x$metadata[ids, n]))) stop("Supplied metadata disagrees with saved ", n, ".")
  }
  md$.cluster <- as.character(md[[cluster_col]])
  observed_clusters <- unique(md$.cluster)
  if (is.null(clusters)) {
    clusters <- if (all(grepl("^[0-9]+$", observed_clusters))) as.character(sort(as.numeric(observed_clusters))) else sort(observed_clusters)
  }
  clusters <- as.character(clusters)
  if (!length(clusters) || anyDuplicated(clusters) || !all(clusters %in% observed_clusters)) stop("Requested clusters must have scored cells and unique labels.")
  md <- md[md$.cluster %in% clusters, , drop = FALSE]
  ids <- rownames(md)
  cv <- x$coverage[x$coverage$status == "scored", , drop = FALSE]
  if (is.null(signatures)) score_names <- rownames(scores) else {
    if (!all(signatures %in% c(cv$signature, cv$score_name))) stop("Unknown requested signature/arm.")
    score_names <- cv$score_name[cv$signature %in% signatures | cv$score_name %in% signatures]
  }
  if (!length(score_names)) stop("No scored signature arms selected.")
  grouping_note <- "Original time labels retained. No T8H-to-T10H display regrouping."
  if (group_by == "condition") {
    if (condition_grouping == "deg_five") {
      md$.condition <- signature_five_conditions(md)
      default_levels <- c("Control", "P4", "N4", "P10", "N10")
      if (!is.null(display_groups)) {
        if (!length(display_groups) || anyDuplicated(display_groups) ||
            !all(display_groups %in% unique(md$.condition))) stop("Display groups must be unique observed approved conditions.")
        default_levels <- default_levels[default_levels %in% display_groups]
        md <- md[md$.condition %in% display_groups, , drop = FALSE]
        ids <- rownames(md)
      }
      if (!is.null(group_levels) && !identical(group_levels, default_levels)) stop("deg_five retains canonical order of explicitly selected conditions.")
      grouping_note <- "Owner-approved DEG-five display: Control/P4/N4/P10/N10; REVO31-N8/T8H negative displayed as N10. Original sample/time/pressure preserved; no score recomputation."
    } else if (is.null(condition_col)) {
      # Literal original metadata: T8H remains T8H, never pooled with T10H.
      md$.condition <- paste(as.character(md$time_point), as.character(md$pressure))
      panel <- unique(md[, c("time_point", "pressure", ".condition")])
      time_number <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(panel$time_point))))
      po <- match(as.character(panel$pressure), c("None", "positive", "negative"))
      default_levels <- unique(panel$.condition[order(time_number, po, panel$.condition, na.last = TRUE)])
    } else {
      if (!condition_col %in% names(md) || anyNA(md[[condition_col]])) stop("Missing condition column.")
      md$.condition <- as.character(md[[condition_col]])
      default_levels <- if (is.factor(md[[condition_col]])) levels(droplevels(md[[condition_col]])) else sort(unique(md$.condition))
      grouping_note <- paste("Custom display grouping from metadata column", condition_col,
        "; exported unchanged with original sample/time/pressure; no internal time recoding.")
    }
    if (is.null(group_levels)) group_levels <- default_levels
    if (anyDuplicated(group_levels) || !all(md$.condition %in% group_levels)) stop("group_levels must include every observed condition exactly once.")
  } else {
    if (is.null(group_levels)) group_levels <- clusters
    if (anyDuplicated(group_levels) || !setequal(group_levels, clusters)) stop("Cluster group_levels must match selected clusters.")
  }
  comparisons <- NULL
  if (!is.null(comparisons_file)) {
    comparisons <- readRDS(comparisons_file)
    cs <- comparisons$settings
    if (score_method != "mean_log_expression" || point_mode != "cells" || !pooled_mean_bar ||
        length(clusters) != 1L || !identical(cs$cluster, clusters) || !identical(cs$cluster_col, cluster_col) ||
        !identical(unname(cs$score_md5), unname(tools::md5sum(scores_file))) ||
        !identical(cs$cell_ids, colnames(scores)) || !identical(cs$method, "paired_t_sample_means") ||
        cs$family_size != 68L || cs$adjustment != "BH" || nrow(comparisons$tests) != 68L ||
        anyDuplicated(paste(comparisons$tests$score_name, comparisons$tests$contrast))) stop("Comparison provenance/family does not match saved scores and requested view.")
    for (g in group_levels) {
      expected_ids <- rownames(x$metadata)[signature_five_conditions(x$metadata) == g]
      if (!setequal(rownames(md)[md$.condition == g], expected_ids)) stop("Comparative views must include ALL scored cells in each selected condition.")
    }
  }
  palette_levels <- if (group_by == "cluster") {
    if (all(grepl("^[0-9]+$", observed_clusters))) as.character(sort(as.numeric(observed_clusters))) else sort(observed_clusters)
  } else if (condition_grouping == "deg_five") {
    c("Control", "P4", "N4", "P10", "N10")
  } else if (is.null(condition_col)) {
    # Existing Control/P4/N4/P10/N10 colours stay fixed; original T8H gets its own colour.
    known <- c("T0H None", "T4H positive", "T4H negative", "T10H positive", "T10H negative")
    c(known, setdiff(group_levels, known))
  } else group_levels
  if (length(palette_levels) > length(viz_palette)) stop("More groups than established palette entries.")
  colors <- setNames(viz_palette[seq_along(palette_levels)], palette_levels)
  if (point_mode == "cells") {
    plotted_groups <- if (group_by == "cluster") md$.cluster else md$.condition
    md$.cell_dot_x <- match(plotted_groups, group_levels) +
      withr::with_seed(point_seed, stats::runif(nrow(md), min = -.06, max = .06))
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  info_path <- file.path(output_root, ".INFO")
  writeLines(c("# Saved signature violin plots", "Status: RUNNING", paste("Started:", Sys.time()),
    paste("Score method:", score_method, "; y-axis:", y_label),
    paste("Scores:", scores_file), paste("Cluster column:", cluster_col), paste("Clusters:", paste(clusters, collapse = ",")),
    paste("View:", group_by, "; group order:", paste(group_levels, collapse = ",")),
    paste("Matched scored cells before cluster selection:", n_matched, "; plotted cells:", nrow(md)),
    paste("Source score cells outside supplied metadata:", sum(!colnames(scores) %in% rownames(metadata))),
    paste("Supplied metadata cells not scored:", sum(!rownames(metadata) %in% colnames(scores))),
    paste("Explicit display groups:", paste(group_levels, collapse = ","), "; selected cell IDs exported; no within-group downsampling."),
    if (is.null(comparisons)) "Scores unchanged; no z-scaling, statistics, activity calls or expression recomputation." else
      paste("Saved comparisons:", comparisons_file, "; paired two-sided t-tests of sample means; global BH family68; no tests recomputed for this view."),
    if (pooled_mean_bar) "Horizontal bar: pooled-cell mean (cell-weighted, descriptive), NOT the equal-pig sample-mean inferential estimand." else "No mean bar.",
    "Violin pools cells descriptively (larger specimens contribute more). No density for singleton or constant groups.",
    if (point_mode == "sample_medians") "White points: sample medians; black point/range: pooled-cell median/IQR." else
      paste("Every dot is one cell; no downsampling, sample-median or pooled-median overlay. Per-plot count equals selected panel cell count. Display-only seed", point_seed,
        "; horizontal jitter Uniform[-0.06,0.06], no vertical jitter; size0.2,alpha0.35,black; ggrastr raster.dpi300,vector text/violins."),
    if (score_method == "mean_log_expression") "Mean = arithmetic average of fixed retained RNA LogNormalize features, zeros included; no [0,1] cap." else "Normalized AUCell AUC; [0,1] bounds.",
    grouping_note), info_path)
  completed <- FALSE
  on.exit(if (!completed) cat("\nStatus: FAILED; incomplete plots.\n", file = info_path, append = TRUE))
  summary_rows <- list()
  sample_rows <- list(); mean_rows <- list(); annotation_rows <- list()
  selected_md <- data.frame(cell_id = rownames(md), md[, unique(c(required, condition_col)), drop = FALSE], check.names = FALSE)
  selected_md$signature_plot_group <- if (group_by == "cluster") md$.cluster else md$.condition
  selected_md$signature_condition_grouping <- if (group_by == "cluster") "cluster" else if (!is.null(condition_col)) paste0("column:", condition_col) else condition_grouping
  if (point_mode == "cells") selected_md$signature_cell_dot_x <- md$.cell_dot_x
  write.csv(selected_md, file.path(output_root, "plotted_cell_metadata.csv"), row.names = FALSE)
  k <- 0L
  for (sn in score_names) {
    coverage <- cv[match(sn, cv$score_name), ]
    all_values <- scores[sn, ids]
    upper <- max(all_values)
    if (upper == 0) upper <- 0.01 # display-only axis span for an all-zero score
    y_limits <- c(0, if (score_method == "AUCell") min(1, upper * 1.08) else upper * 1.08)
    panels <- if (group_by == "cluster") list(all = ids) else setNames(lapply(clusters,
      function(cl) rownames(md)[md$.cluster == cl]), clusters)
    for (panel_name in names(panels)) {
      cc <- panels[[panel_name]]
      d <- data.frame(cell_id = cc, value = as.numeric(scores[sn, cc]),
        sample = as.character(md[cc, "sample"]), pig = as.character(md[cc, "pig"]),
        group = if (group_by == "cluster") md[cc, ".cluster"] else md[cc, ".condition"], stringsAsFactors = FALSE)
      d$group <- factor(d$group, levels = group_levels)
      if (point_mode == "cells") d$cell_dot_x <- md[cc, ".cell_dot_x"]
      present_groups <- group_levels[group_levels %in% as.character(d$group)]
      summaries <- do.call(rbind, lapply(present_groups, function(g) {
        z <- d[d$group == g, ]
        q <- stats::quantile(z$value, c(.25, .5, .75), names = FALSE)
        data.frame(group = g, n_cells = nrow(z), n_samples = length(unique(z$sample)),
          n_pigs = length(unique(z$pig)), q25 = q[1], median = q[2], q75 = q[3],
          draw_density = nrow(z) >= 2 && length(unique(z$value)) >= 2)
      }))
      means <- stats::aggregate(value ~ group, data = d, FUN = mean)
      means$x <- match(as.character(means$group), group_levels)
      ann <- NULL
      if (!is.null(comparisons)) {
        tt <- comparisons$tests
        ann <- tt[tt$score_name == sn & tt$first %in% group_levels & tt$second %in% group_levels, , drop = FALSE]
        if (!nrow(ann) || any(!is.finite(ann$p_adj))) stop("No valid predeclared comparisons for selected display groups.")
        ann$x <- match(ann$first, group_levels); ann$xend <- match(ann$second, group_levels)
        ann <- ann[order(abs(ann$x - ann$xend), ann$contrast), , drop = FALSE]
        ann$y <- upper * (1.08 + .15 * seq_len(nrow(ann)))
        ann$label <- paste0("BH-adjusted p = ", format.pval(ann$p_adj, digits = 3, eps = 0))
        y_limits <- c(0, max(ann$y) + upper * .16)
      }
      samples <- stats::aggregate(value ~ group + sample + pig, data = d, FUN = stats::median)
      # Deterministic small offsets expose sample points without consuming scoring RNG.
      sample_order <- sort(unique(as.character(md$sample)), method = "radix")
      offsets <- setNames(seq(-.16, .16, length.out = length(sample_order)), sample_order)
      samples$x <- match(as.character(samples$group), group_levels) + offsets[samples$sample]
      summaries$group <- factor(summaries$group, levels = group_levels)
      density_groups <- as.character(summaries$group[summaries$draw_density])
      title <- gsub("_", " ", coverage$signature)
      subtitle <- paste0(coverage$arm, " arm | ", coverage$effective_pig_genes, " unique pig features; ",
        coverage$human_in_universe, "/", coverage$human_total, " human members represented",
        if (group_by == "condition") paste0(" | cluster ", panel_name) else "")
      labels <- setNames(paste0(summaries$group, "\n", summaries$n_cells, " cells",
        if (point_mode == "sample_medians") paste0("\n", summaries$n_samples, " samples\n", summaries$n_pigs, " pigs") else ""), as.character(summaries$group))
      # Retain absent conditions on the axis; do not fabricate density/zero observations.
      missing_groups <- setdiff(group_levels, names(labels))
      labels <- c(labels, setNames(paste0(missing_groups, "\n0 cells"), missing_groups))
      p <- ggplot2::ggplot(d, ggplot2::aes(x = group, y = value, fill = group)) +
        (if (length(density_groups)) ggplot2::geom_violin(data = d[as.character(d$group) %in% density_groups, , drop = FALSE],
          scale = "width", trim = TRUE, width = .8, linewidth = .25, alpha = .85)) +
        (if (point_mode == "sample_medians") ggplot2::geom_linerange(data = summaries, ggplot2::aes(x = group, ymin = q25, ymax = q75),
          inherit.aes = FALSE, linewidth = .7, color = "black")) +
        (if (point_mode == "sample_medians") ggplot2::geom_point(data = summaries, ggplot2::aes(x = group, y = median),
          inherit.aes = FALSE, size = 1.8, color = "black")) +
        (if (point_mode == "sample_medians") ggplot2::geom_point(data = samples, ggplot2::aes(x = x, y = value),
          inherit.aes = FALSE, shape = 21, fill = "white", color = "black", stroke = .3, size = 1.6)) +
        (if (point_mode == "cells") ggrastr::geom_point_rast(data = d, ggplot2::aes(x = cell_dot_x, y = value),
          inherit.aes = FALSE, color = "black", size = .2, alpha = .35, raster.dpi = 300)) +
        (if (pooled_mean_bar) ggplot2::geom_segment(data = means, ggplot2::aes(x = x - .18, xend = x + .18, y = value, yend = value), inherit.aes = FALSE, color = "white", linewidth = 2)) +
        (if (pooled_mean_bar) ggplot2::geom_segment(data = means, ggplot2::aes(x = x - .18, xend = x + .18, y = value, yend = value), inherit.aes = FALSE, color = "black", linewidth = .8)) +
        (if (!is.null(ann)) ggplot2::geom_segment(data = ann, ggplot2::aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = .3)) +
        (if (!is.null(ann)) ggplot2::geom_segment(data = ann, ggplot2::aes(x = x, xend = x, y = y, yend = y - upper * .025), inherit.aes = FALSE, linewidth = .3)) +
        (if (!is.null(ann)) ggplot2::geom_segment(data = ann, ggplot2::aes(x = xend, xend = xend, y = y, yend = y - upper * .025), inherit.aes = FALSE, linewidth = .3)) +
        (if (!is.null(ann)) ggplot2::geom_text(data = ann, ggplot2::aes(x = (x + xend)/2, y = y + upper * .025, label = label), inherit.aes = FALSE, size = 2.7, vjust = 0)) +
        (if (length(density_groups)) ggplot2::scale_fill_manual(values = colors, drop = FALSE)) +
        ggplot2::scale_x_discrete(limits = group_levels, labels = labels[group_levels], drop = FALSE) +
        ggplot2::scale_y_continuous(limits = y_limits, expand = ggplot2::expansion(mult = c(.02, .03))) +
        ggplot2::labs(title = deg_wrap(title, 70), subtitle = deg_wrap(subtitle, 95),
          x = if (group_by == "cluster") "Cluster" else if (condition_grouping == "deg_five") "Condition" else if (!is.null(condition_col)) condition_col else "Original time / pressure condition",
          y = y_label,
          caption = deg_wrap(paste0(if (!is.null(ann)) paste0(
            "Dots: all cells; bar: pooled-cell mean (cell-weighted display only). Tests: two-sided paired t-tests of sample means, equal weight per pig; n = ",
            paste(unique(ann$n_pairs), collapse = "/"), " matched pigs. BH adjusted jointly across all 68 tests; bars are not the test estimand.") else if (point_mode == "cells")
            paste0("Each black dot is one cell; no downsampling. ", if (pooled_mean_bar) "Bar: pooled-cell mean (cell-weighted display). " else "",
              "Descriptive pooled distributions; cells are not independent replicates. No statistical tests.") else
            "Descriptive cell distributions; black: median/IQR; white points: sample medians. Larger samples contribute more cells. No statistical tests.",
            if (score_method == "mean_log_expression") " Score: mean of retained RNA LogNormalize values, including zeros." else "",
            if (condition_grouping == "deg_five" && "N10" %in% group_levels) " N10 includes REVO31-N8 (original T8H)." else ""), 100)) +
        deg_theme() + ggplot2::theme(legend.position = "none",
          axis.text.x = ggplot2::element_text(size = 8, angle = if (point_mode == "cells") 45 else 0,
            hjust = if (point_mode == "cells") 1 else .5))
      stem <- paste0(sn, if (group_by == "condition") paste0("__cluster_", panel_name) else "__clusters")
      stem <- gsub("[^A-Za-z0-9_.-]", "_", stem)
      save_deg_figure(p, output_root, stem, width = max(170, length(group_levels) * 26), height = if (is.null(ann)) 135 else 155)
      k <- k + 1L
      summaries$score_name <- sn; summaries$panel <- panel_name
      samples$score_name <- sn; samples$panel <- panel_name
      summary_rows[[k]] <- summaries; sample_rows[[k]] <- samples
      if (pooled_mean_bar) { means$score_name <- sn; means$panel <- panel_name; mean_rows[[k]] <- means }
      if (!is.null(ann)) { ann$panel <- panel_name; annotation_rows[[k]] <- ann }
    }
  }
  mean_bars <- if (length(mean_rows)) do.call(rbind, mean_rows) else NULL
  annotations <- if (length(annotation_rows)) do.call(rbind, annotation_rows) else NULL
  if (!is.null(mean_bars)) { names(mean_bars)[names(mean_bars) == "value"] <- "pooled_cell_mean"; write.csv(mean_bars, file.path(output_root, "pooled_cell_mean_bars.csv"), row.names = FALSE) }
  if (!is.null(annotations)) write.csv(annotations, file.path(output_root, "comparison_annotations.csv"), row.names = FALSE)
  summaries <- do.call(rbind, summary_rows)
  samples <- do.call(rbind, sample_rows)
  write.csv(summaries, file.path(output_root, "group_summaries.csv"), row.names = FALSE)
  write.csv(samples, file.path(output_root, "sample_medians.csv"), row.names = FALSE)
  stopifnot(nrow(read.csv(file.path(output_root, "sample_medians.csv"))) == nrow(samples))
  cat("\nStatus: COMPLETED; summary tables reopened. Saved figures require visual review.\n", file = info_path, append = TRUE)
  completed <- TRUE
  invisible(list(group_summaries = summaries, sample_medians = samples, cell_ids = ids,
    score_method = score_method, point_mode = point_mode, mean_bars = mean_bars, annotations = annotations,
    cell_dot_ids = if (point_mode == "cells") ids else character()))
}
