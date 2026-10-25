# Approved 17-arm, four-contrast family on SAVED mean log-expression scores.
# Biological unit: one sample mean per pig/condition; never individual cells.
# source("R/DEG/compare_signature_scores.R")
# compare_signature_scores("<mean score root>/signature_scores.rds", "<NEW comparison directory>")
source("R/DEG/signature_conditions.R")

compare_signature_scores <- function(scores_file, output_root,
    cluster_col = "harmony_snn_res.0.3", cluster = "1") {
  if (file.exists(output_root)) stop("output_root exists; use a new directory.")
  contrasts <- data.frame(contrast = c("P4-N4", "P10-N10", "P4-Control", "N4-Control"),
    first = c("P4", "P10", "P4", "N4"), second = c("N4", "N10", "Control", "Control"),
    expected_pairs = c(5L, 4L, 5L, 5L), stringsAsFactors = FALSE)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  info <- file.path(output_root, ".INFO")
  writeLines(c("# Paired sample-mean signature comparisons", "Status: RUNNING", paste("Started:", Sys.time()),
    paste("Read-only saved scores:", scores_file), paste("Cluster:", cluster_col, cluster),
    "Predeclared family:17arms x4contrasts=P4-N4,P10-N10,P4-Control,N4-Control; exactly68tests.",
    "Two-sided paired t-tests on arithmetic per-sample means, matched by pig; equal pig weights, not cell tests.",
    "Global BH across all68; FDR<0.05. Ordinary paired95%CI is NOT multiplicity-adjusted.",
    "One sample per pig-condition; require complete5/4/5/5pairs; no pair deletion, outlier filtering or normality gates.",
    "REVO31-N8 originalT8H/negative displayedN10; preserved source identities. Separate arms, no composite.",
    "Degenerate paired differences stop the entire family; no fabricated p-values or reduced-family correction."), info)
  completed <- FALSE
  on.exit(if (!completed) cat("\nStatus: FAILED; no validated complete68test family.\n", file = info, append = TRUE))
  score_md5 <- tools::md5sum(scores_file)
  x <- readRDS(scores_file); md <- x$metadata; scores <- x$scores
  required <- c(cluster_col, "sample", "pig", "time_point", "pressure")
  if (!identical(x$settings$method, "mean_log_expression") || !all(required %in% names(md)) ||
      anyNA(md[, required, drop = FALSE]) || anyDuplicated(rownames(md)) ||
      !identical(colnames(scores), rownames(md)) || nrow(scores) != 17L || anyDuplicated(rownames(scores)) ||
      any(!is.finite(scores)) || any(scores < 0) || length(cluster) != 1L ||
      !all(as.character(md[[cluster_col]]) == cluster) || !identical(x$settings$cluster_col, cluster_col) ||
      !identical(x$settings$clusters, cluster)) stop("Require complete saved mean scores for exactly the requested source cluster and17arms.")
  coverage <- x$coverage[match(rownames(scores), x$coverage$score_name), , drop = FALSE]
  if (anyNA(coverage$score_name) || any(coverage$status != "scored")) stop("Missing original arm coverage.")
  md$condition <- signature_five_conditions(md)
  panel <- unique(md[, c("sample", "pig", "time_point", "pressure", "condition"), drop = FALSE])
  panel <- panel[order(as.character(panel$sample), method = "radix"), , drop = FALSE]
  if (nrow(panel) != 23L || anyDuplicated(panel$sample) ||
      anyDuplicated(paste(panel$pig, panel$condition, sep = "\t"))) stop("Require original23samples, one metadata tuple/sample and one sample/pig-condition.")
  pairs <- lapply(seq_len(nrow(contrasts)), function(i) {
    a <- panel[panel$condition == contrasts$first[i], ]; b <- panel[panel$condition == contrasts$second[i], ]
    pigs <- sort(as.character(a$pig), method = "radix")
    if (!setequal(pigs, as.character(b$pig)) || length(pigs) != contrasts$expected_pairs[i]) stop("Incomplete approved pig pairing for ", contrasts$contrast[i], "; do not silently drop pigs.")
    data.frame(contrast = contrasts$contrast[i], pig = pigs,
      first_sample = as.character(a$sample[match(pigs, a$pig)]),
      second_sample = as.character(b$sample[match(pigs, b$pig)]), stringsAsFactors = FALSE)
  })
  pairs <- do.call(rbind, pairs)
  panel$n_cells <- vapply(panel$sample, function(s) sum(as.character(md$sample) == s), integer(1))
  sample_values <- vapply(as.character(panel$sample), function(s)
    rowMeans(scores[, as.character(md$sample) == s, drop = FALSE]), numeric(nrow(scores)))
  rownames(sample_values) <- rownames(scores)
  sample_means <- do.call(rbind, lapply(rownames(scores), function(sn)
    data.frame(score_name = sn, panel, mean_score = as.numeric(sample_values[sn, ]), row.names = NULL)))
  condition_means <- stats::aggregate(mean_score ~ score_name + condition, sample_means, mean)
  names(condition_means)[3] <- "mean_of_sample_means"
  results <- list(); k <- 0L
  for (sn in rownames(scores)) for (i in seq_len(nrow(contrasts))) {
    pair <- pairs[pairs$contrast == contrasts$contrast[i], ]
    a <- as.numeric(sample_values[sn, pair$first_sample]); b <- as.numeric(sample_values[sn, pair$second_sample])
    differences <- a - b
    if (!is.finite(stats::sd(differences)) || stats::sd(differences) == 0) stop("Degenerate paired differences: ", sn, " ", contrasts$contrast[i])
    tt <- tryCatch(stats::t.test(a, b, paired = TRUE, alternative = "two.sided", conf.level = .95),
      error = function(e) stop(sn, " ", contrasts$contrast[i], ": ", conditionMessage(e)))
    k <- k + 1L
    results[[k]] <- data.frame(score_name = sn, contrast = contrasts$contrast[i],
      first = contrasts$first[i], second = contrasts$second[i], n_pairs = nrow(pair),
      mean_first = mean(a), mean_second = mean(b), paired_mean_difference = mean(differences),
      ci_low = tt$conf.int[1], ci_high = tt$conf.int[2], t = unname(tt$statistic),
      df = unname(tt$parameter), p_value = tt$p.value, stringsAsFactors = FALSE)
  }
  tests <- do.call(rbind, results)
  stopifnot(nrow(tests) == 68L, all(is.finite(as.matrix(tests[, c("paired_mean_difference", "ci_low", "ci_high", "t", "df", "p_value")]))))
  tests$p_adj <- stats::p.adjust(tests$p_value, method = "BH", n = 68L)
  tests$significant_fdr05 <- tests$p_adj < .05
  tests$effective_pig_genes <- coverage$effective_pig_genes[match(tests$score_name, coverage$score_name)]
  tests$small_set_warning <- coverage$small_set_warning[match(tests$score_name, coverage$score_name)]
  result <- list(tests = tests, sample_means = sample_means, condition_means = condition_means,
    pairs = pairs, coverage = coverage,
    settings = list(scores_file = scores_file, score_md5 = score_md5, cell_ids = colnames(scores),
      cluster_col = cluster_col, cluster = cluster, method = "paired_t_sample_means", alternative = "two.sided",
      contrasts = contrasts, family_size = 68L, adjustment = "BH", fdr = .05, conf_level = .95,
      ci_multiplicity_adjusted = FALSE, grouping = "deg_five", session = utils::sessionInfo()))
  for (n in c("tests", "sample_means", "condition_means", "pairs")) write.csv(result[[n]], file.path(output_root, paste0(n, ".csv")), row.names = FALSE)
  saveRDS(result, file.path(output_root, "comparisons.rds"))
  stopifnot(identical(result, readRDS(file.path(output_root, "comparisons.rds"))), identical(score_md5, tools::md5sum(scores_file)))
  cat("\nStatus: COMPLETED;68finite tests; global BH; exact RDS reopen and score-file hash unchanged.\n",
    "Plots may only reuse this family; repeated brackets do not create new tests.\n", file = info, append = TRUE)
  completed <- TRUE
  invisible(result)
}
