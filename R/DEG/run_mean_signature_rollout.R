# Fixed authorized rollout orchestration; scientific methods remain in existing helpers.
# Source from project root, then run_mean_signature_rollout(<preflight directory>, <NEW run-record directory>).
# Sequential cluster execution; original Endothelial C1 is never regenerated.
source("R/DEG/run_mean_signature_scores.R")
source("R/DEG/compare_signature_scores.R")
source("R/plots/plot_signature_scores.R")
source("R/DEG/preflight_mean_signature_rollout.R")

run_mean_signature_rollout <- function(preflight_root, run_root) {
  if (file.exists(run_root)) stop("Use a NEW run-record directory.")
  pf <- readRDS(file.path(preflight_root, "preflight.rds"))
  stopifnot(identical(pf$manifest, mean_signature_rollout_manifest()), all(pf$pairs$full_original_pairing))
  protected <- readRDS(file.path(preflight_root, "protected_md5.rds"))
  stopifnot(identical(protected, tools::md5sum(names(protected))))
  manifest <- pf$manifest; todo <- manifest[manifest$action != "preserve_completed", ]
  stopifnot(nrow(todo) == 27L, !any(file.exists(todo$output_root)))
  dir.create(run_root, recursive = TRUE)
  run_info <- file.path(run_root, ".INFO")
  writeLines(c("# Authorized mean-signature rollout", "Status: RUNNING",
    "27new clusters; completed EndothelialC1 preserved. Sequential existing methods, per-cluster68BH only.",
    "Runtime degenerate inference remains pending owner decision; no fabricated tests. PDF/PNG render review follows numerical execution."), run_info)
  write.csv(manifest, file.path(run_root, "manifest.csv"), row.names = FALSE)
  reference <- readRDS(file.path(manifest$output_root[manifest$action == "preserve_completed"], "signature_scores.rds"))
  mapping <- "results/misc/Ensembl_ID/signature-scoring-v2/pig_mapping.csv"
  universe <- "results/misc/Ensembl_ID/signature-scoring/common_count_universe.rds"
  original_save <- save_deg_figure
  on.exit(assign("save_deg_figure", original_save, envir = .GlobalEnv), add = TRUE)
  inventory <- list(); figure_checks <- list()
  groups <- c("Control", "P4", "N4", "P10", "N10")
  layouts <- list(five_conditions = groups, P4_N4 = c("P4", "N4"), P10_N10 = c("P10", "N10"), Control_P4_N4 = c("Control", "P4", "N4"))
  for (population in unique(todo$population)) {
    units <- todo[todo$population == population, ]
    for (i in seq_len(nrow(units))) {
      unit <- units[i, ]; key <- paste0(population, "_cluster_", unit$cluster)
      log <- file(file.path(run_root, paste0(key, ".log")), open = "wt")
      sink(log); sink(log, type = "message")
      outcome <- tryCatch({
        started <- Sys.time()
        a <- run_mean_signature_scores(unit$input, mapping, universe, unit$output_root,
          cluster_col = unit$cluster_col, clusters = unit$cluster)
        gc()
        # Independent original-RDS arithmetic/identity validation, all scored cells.
        object <- readRDS(unit$input); md0 <- object[[]]
        expected_ids <- rownames(md0)[as.character(md0[[unit$cluster_col]]) == unit$cluster]
        expected_row <- pf$clusters[pf$clusters$population == population & pf$clusters$cluster == unit$cluster, ]
        stopifnot(identical(dim(a$scores), c(17L, as.integer(expected_row$n_cells))),
          setequal(colnames(a$scores), expected_ids), !anyDuplicated(colnames(a$scores)),
          identical(a$metadata, md0[colnames(a$scores), , drop = FALSE]),
          identical(a$gene_sets, reference$gene_sets), identical(a$mapping, reference$mapping), identical(a$coverage, reference$coverage),
          identical(a, readRDS(file.path(unit$output_root, "signature_scores.rds"))))
        max_score_error <- 0; checked_values <- 0L
        for (layer in SeuratObject::Layers(object[["RNA"]], search = "^data[.]")) {
          sample <- sub("^data[.]", "", layer)
          ids <- colnames(a$scores)[as.character(a$metadata$sample) == sample]
          mat <- SeuratObject::LayerData(object[["RNA"]], layer = layer, features = unique(unlist(a$gene_sets)), cells = ids)
          for (sn in names(a$gene_sets)) {
            manual <- colSums(as.matrix(mat[a$gene_sets[[sn]], ids, drop = FALSE])) / length(a$gene_sets[[sn]])
            err <- max(abs(manual - a$scores[sn, ids])); stopifnot(err < 1e-12)
            max_score_error <- max(max_score_error, err); checked_values <- checked_values + length(ids)
          }
        }
        stopifnot(checked_values == length(a$scores))
        rm(object, md0, mat); gc()
        md <- a$metadata
        cond <- ifelse(md$time_point == "T0H", "Control", ifelse(md$time_point == "T4H",
          ifelse(md$pressure == "positive", "P4", "N4"), ifelse(md$pressure == "positive", "P10", "N10")))
        expected_samples <- pf$samples[pf$samples$population == population & pf$samples$cluster == unit$cluster, ]
        for (j in seq_len(nrow(expected_samples))) stopifnot(sum(md$sample == expected_samples$sample[j]) == expected_samples$n_cells[j])
        stopifnot(identical(cond, signature_five_conditions(md)), all(md$time_point[md$sample == "REVO31-N8"] == "T8H"))
        sf <- file.path(unit$output_root, "signature_scores.rds")
        comparison_root <- file.path(unit$output_root, "comparisons_paired_sample_means")
        z <- tryCatch(compare_signature_scores(sf, comparison_root, unit$cluster_col, unit$cluster), error = function(e) {
          if (!grepl("Degenerate paired differences|data are essentially constant", conditionMessage(e))) stop(e)
          structure(list(error = conditionMessage(e)), class = "inference_pending")
        })
        pending <- inherits(z, "inference_pending")
        max_stat_error <- NA_real_; max_bh_error <- NA_real_
        if (!pending) {
          sm <- z$sample_means
          for (j in seq_len(nrow(sm))) {
            r <- sm[j, ]; ids <- rownames(md)[as.character(md$sample) == r$sample]
            stopifnot(length(ids) == r$n_cells, abs(sum(a$scores[r$score_name, ids])/length(ids) - r$mean_score) < 1e-14)
          }
          max_stat_error <- 0
          for (j in seq_len(nrow(z$tests))) {
            r <- z$tests[j, ]; aa <- sm[sm$score_name == r$score_name & sm$condition == r$first, ]; bb <- sm[sm$score_name == r$score_name & sm$condition == r$second, ]
            stopifnot(setequal(aa$pig, bb$pig), nrow(aa) == r$n_pairs)
            d <- aa$mean_score - bb$mean_score[match(aa$pig, bb$pig)]; n <- length(d); delta <- sum(d)/n
            se <- sqrt(sum((d-delta)^2)/(n-1)/n); t <- delta/se
            ci <- delta + c(-1,1)*stats::qt(.975,n-1)*se; p <- 2*stats::pt(-abs(t),n-1)
            err <- max(abs(c(delta,t,p,ci) - c(r$paired_mean_difference,r$t,r$p_value,r$ci_low,r$ci_high)))
            stopifnot(err < 1e-10, r$df == n-1); max_stat_error <- max(max_stat_error, err)
          }
          p <- z$tests$p_value; o <- order(p); bh <- numeric(68)
          bh[o] <- pmin(1, rev(cummin(rev(p[o]*68/seq_along(p)))))
          max_bh_error <- max(abs(bh-z$tests$p_adj))
          stopifnot(nrow(z$tests) == 68L, max_bh_error < 1e-14,
            identical(z$tests$significant_fdr05, z$tests$p_adj < .05), identical(z, readRDS(file.path(comparison_root, "comparisons.rds"))))
          for (nm in c("tests", "sample_means", "condition_means", "pairs")) {
            disk <- read.csv(file.path(comparison_root, paste0(nm, ".csv")))
            stopifnot(nrow(disk) == nrow(z[[nm]]), identical(names(disk), names(z[[nm]])))
            for (cn in names(disk)) if (is.numeric(disk[[cn]])) stopifnot(isTRUE(all.equal(as.numeric(disk[[cn]]), as.numeric(z[[nm]][[cn]]), tolerance = 1e-12))) else stopifnot(identical(as.character(disk[[cn]]), as.character(z[[nm]][[cn]])))
          }
        } else {
          # Report only, no invented p-values or altered family. Independent units continue.
          cat("INFERENCE_PENDING: ", z$error, "\n", file = file.path(unit$output_root, ".INFO"), append = TRUE)
          writeLines(c("# Inference pending owner decision", z$error,
            "Scores validated. No comparison figures/fabricated p-values/reduced-family BH. Five-condition descriptive view only."), file.path(comparison_root, "PENDING.INFO"))
        }
        protected_results <- c(sf, if (!pending) file.path(comparison_root, c("comparisons.rds", "tests.csv", "sample_means.csv", "condition_means.csv", "pairs.csv")))
        result_hashes <- tools::md5sum(protected_results)
        current_layout <- NULL
        # Validate every built figure at the existing save boundary; no method changes.
        assign("save_deg_figure", function(figure, directory, stem, width = 170, height) {
          sn <- sub("__cluster_[^_]+$", "", stem)
          point <- which(vapply(figure$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))
          stopifnot(length(point) == 1L)
          d <- figure$layers[[point]]$data; selected_groups <- levels(d$group)
          ids <- rownames(md)[cond %in% selected_groups]; built <- ggplot2::ggplot_build(figure)
          stopifnot(setequal(d$cell_id, ids), !anyDuplicated(d$cell_id), identical(d$value, as.numeric(a$scores[sn,d$cell_id])),
            nrow(built$data[[point]]) == length(ids), !anyNA(built$data[[point]]$y), identical(built$data[[point]]$y, d$value),
            max(d$value) < figure$scales$get_scales("y")$limits[2])
          n_bars <- 0L; n_brackets <- 0L
          if (current_layout != "five_conditions") {
            ti <- which(vapply(figure$layers, function(l) inherits(l$geom, "GeomText"), logical(1))); stopifnot(length(ti) == 1L)
            ann <- figure$layers[[ti]]$data; tt <- z$tests[z$tests$score_name == sn & z$tests$first %in% selected_groups & z$tests$second %in% selected_groups, ]
            stopifnot(setequal(ann$contrast, tt$contrast)); tt <- tt[match(ann$contrast,tt$contrast), ]
            stopifnot(identical(ann$p_adj,tt$p_adj), identical(ann$x,match(tt$first,selected_groups)), identical(ann$xend,match(tt$second,selected_groups)),
              identical(built$data[[ti]]$label, paste0("BH-adjusted p = ",format.pval(tt$p_adj,digits=3,eps=0))),
              max(built$data[[ti]]$y) < figure$scales$get_scales("y")$limits[2], max(d$value) < min(ann$y))
            segments <- which(vapply(figure$layers,function(l) inherits(l$geom,"GeomSegment"),logical(1))); stopifnot(length(segments) == 5L)
            bars <- figure$layers[[segments[1]]]$data
            for (j in seq_len(nrow(bars))) {
              ids <- rownames(md)[cond == as.character(bars$group[j])]
              sm <- z$sample_means[z$sample_means$score_name == sn & z$sample_means$condition == as.character(bars$group[j]), ]
              stopifnot(identical(bars$value[j],mean(a$scores[sn,ids])), abs(bars$value[j]-sum(sm$mean_score*sm$n_cells)/sum(sm$n_cells)) < 1e-14)
            }
            stopifnot(identical(built$data[[segments[1]]]$y,bars$value), identical(built$data[[segments[2]]]$y,bars$value),
              all(built$data[[segments[3]]]$x == ann$x), all(built$data[[segments[3]]]$xend == ann$xend), all(built$data[[segments[3]]]$y == ann$y),
              all(built$data[[segments[4]]]$x == ann$x), all(built$data[[segments[5]]]$x == ann$xend))
            n_bars <- nrow(bars); n_brackets <- nrow(ann)
          } else stopifnot(all(vapply(figure$layers, function(l)
            inherits(l$geom, "GeomPoint") || inherits(l$geom, "GeomViolin"), logical(1))))
          figure_checks[[length(figure_checks)+1L]] <<- data.frame(population=population,cluster=unit$cluster,layout=current_layout,score_name=sn,
            n_cell_dots=nrow(d),n_mean_bars=n_bars,n_brackets=n_brackets,pdf=file.path(directory,paste0(stem,".pdf")),png=file.path(directory,paste0(stem,".png")))
          original_save(figure,directory,stem,width=width,height=height)
        }, envir = .GlobalEnv)
        selected_layouts <- if (pending) "five_conditions" else names(layouts)
        for (layout in selected_layouts) {
          current_layout <- layout; descriptive <- layout == "five_conditions"
          plot_root <- if (descriptive) file.path(unit$output_root,"violins_conditions") else file.path(comparison_root,layout)
          p <- plot_signature_scores(sf,unit$cluster_col,plot_root,clusters=unit$cluster,group_by="condition",condition_grouping="deg_five",point_mode="cells",
            display_groups=if (descriptive) NULL else layouts[[layout]],comparisons_file=if (descriptive) NULL else file.path(comparison_root,"comparisons.rds"),pooled_mean_bar=!descriptive)
          expected <- rownames(md)[cond %in% layouts[[layout]]]
          exported <- read.csv(file.path(plot_root,"plotted_cell_metadata.csv"),check.names=FALSE)
          stopifnot(setequal(p$cell_dot_ids,expected),identical(exported$cell_id,p$cell_dot_ids))
          for (cn in c("sample","pig","time_point","pressure",unit$cluster_col)) stopifnot(identical(as.character(exported[[cn]]),as.character(md[exported$cell_id,cn])))
          offset <- exported$signature_cell_dot_x - match(exported$signature_plot_group,layouts[[layout]])
          stopifnot(isTRUE(all.equal(offset,withr::with_seed(1234L,runif(nrow(exported),-.06,.06)),tolerance=1e-12)))
          for (nm in c("group_summaries","sample_medians",if (!descriptive) c("pooled_cell_mean_bars","comparison_annotations"))) {
            disk <- read.csv(file.path(plot_root,paste0(nm,".csv")))
            memory <- switch(nm,group_summaries=p$group_summaries,sample_medians=p$sample_medians,pooled_cell_mean_bars=p$mean_bars,comparison_annotations=p$annotations)
            stopifnot(nrow(disk)==nrow(memory),identical(names(disk),names(memory)))
            for (cn in names(disk)) if (is.numeric(disk[[cn]])) stopifnot(isTRUE(all.equal(as.numeric(disk[[cn]]),as.numeric(memory[[cn]]),tolerance=1e-12))) else stopifnot(identical(as.character(disk[[cn]]),as.character(memory[[cn]])))
          }
          stopifnot(length(list.files(plot_root,pattern="[.]pdf$"))==17L,length(list.files(plot_root,pattern="[.]png$"))==17L)
        }
        assign("save_deg_figure",original_save,envir=.GlobalEnv)
        stopifnot(identical(result_hashes,tools::md5sum(protected_results)))
        cat(paste0("\nRollout numerical/built-plot validation PASS: ",checked_values,"score values; max arithmetic error",max_score_error,
          "; stats error",max_stat_error,"; BH error",max_bh_error,". Every selected cell/bar/bracket verified; saved scores/tests unchanged through plotting.\n",
          "Figure reopen/render review pending. Per-cluster68BH only; no study-wide FDR guarantee. Sparse samples retained without filtering.\n"),file=file.path(unit$output_root,".INFO"),append=TRUE)
        data.frame(population=population,cluster=unit$cluster,n_cells=ncol(a$scores),n_samples=length(unique(md$sample)),
          score_values_checked=checked_values,max_score_error=max_score_error,max_stat_error=max_stat_error,max_bh_error=max_bh_error,
          n_tests=if (pending) NA_integer_ else 68L,n_pdf=17L*length(selected_layouts),n_png=17L*length(selected_layouts),
          sig_P4_N4=if (pending) NA_integer_ else sum(z$tests$significant_fdr05[z$tests$contrast=="P4-N4"]),
          sig_P10_N10=if (pending) NA_integer_ else sum(z$tests$significant_fdr05[z$tests$contrast=="P10-N10"]),
          sig_P4_Control=if (pending) NA_integer_ else sum(z$tests$significant_fdr05[z$tests$contrast=="P4-Control"]),
          sig_N4_Control=if (pending) NA_integer_ else sum(z$tests$significant_fdr05[z$tests$contrast=="N4-Control"]),
          output_root=unit$output_root,status=if (pending) "INFERENCE_PENDING" else "NUMERICAL_PLOT_RECORDS_PASS_RENDER_PENDING",
          error=if (pending) z$error else "",elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")))
      }, error = function(e) e)
      sink(type="message"); sink(); close(log)
      assign("save_deg_figure",original_save,envir=.GlobalEnv)
      if (inherits(outcome,"error")) {
        writeLines(conditionMessage(outcome),file.path(run_root,paste0(key,"_FAILED.INFO")))
        cat("\nStatus: FAILED; ",key,": ",conditionMessage(outcome),"\n",file=run_info,append=TRUE)
        stop("Unexpected failure in ",key,": ",conditionMessage(outcome))
      }
      inventory[[length(inventory)+1L]] <- outcome
      write.csv(do.call(rbind,inventory),file.path(run_root,"cluster_inventory.csv"),row.names=FALSE)
      write.csv(do.call(rbind,figure_checks),file.path(run_root,"figure_inventory.csv"),row.names=FALSE)
      if (outcome$status == "INFERENCE_PENDING") message("INFERENCE_PENDING ",key,": ",outcome$error)
      rm(a,z,md); gc()
    }
    stopifnot(identical(protected,tools::md5sum(names(protected))))
    message("SUBSET_NUMERICAL_COMPLETE ",population," units=",nrow(units),"; original sources/EndothelialC1 unchanged; render review pending.")
    cat("\nSubset numerical phase completed: ",population,"; units=",nrow(units),"; see inventory for pending inference. Original protected hashes unchanged.\n",file=run_info,append=TRUE)
  }
  result <- list(clusters=do.call(rbind,inventory),figures=do.call(rbind,figure_checks))
  saveRDS(result,file.path(run_root,"rollout_inventory.rds"))
  stopifnot(nrow(result$clusters)==27L,sum(result$clusters$n_cells)==111775L,
    identical(result,readRDS(file.path(run_root,"rollout_inventory.rds"))),
    identical(protected,tools::md5sum(names(protected))))
  cat("\nStatus: NUMERICAL_PHASE_COMPLETED; exact inventory reopen/protected hashes PASS. Check inventory for inference_pending; all PDF/PNG render review and final tech-oracle inspection still required.\n",file=run_info,append=TRUE)
  invisible(result)
}
