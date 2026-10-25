# Shared frozen pig-mapping preparation for AUCell and mean log-expression.
# No expression values, scoring methods or eligibility thresholds are chosen here.
prepare_signature_gene_sets <- function(mapping_file, features, signatures) {
  mapping <- read.csv(mapping_file, stringsAsFactors = FALSE, check.names = FALSE)
  membership <- do.call(rbind, lapply(names(signatures), function(n) {
    do.call(rbind, lapply(c("positive", "negative"), function(a) {
      if (length(signatures[[n]][[a]])) data.frame(signature = n, arm = a,
        human_symbol = signatures[[n]][[a]], stringsAsFactors = FALSE)
    }))
  }))
  key <- function(d) paste(d$signature, d$arm, d$human_symbol, sep = "\t")
  if (anyDuplicated(key(mapping)) || !setequal(key(mapping), key(membership))) {
    stop("Mapping must account for every original signed membership exactly once.")
  }
  mapping <- mapping[match(key(membership), key(mapping)), , drop = FALSE]
  if (anyNA(mapping$approved) || anyNA(mapping$feature) || anyNA(mapping$reason)) stop("Incomplete mapping decision.")
  retained <- mapping$approved
  # Earlier audit flags include rejected alternatives as well as genuinely unresolved
  # identities. Do not discard a resolved stable target solely for its stale alias.
  unresolved <- if ("unresolved_identity" %in% names(mapping)) mapping$unresolved_identity else nzchar(mapping$conflicted_features)
  if (!is.logical(retained) || !is.logical(unresolved) || anyNA(unresolved) ||
      any(retained & (unresolved | grepl(";", mapping$target_ids) | grepl(";", mapping$supported_features)))) {
    stop("Approved mapping violates conservative single-target/single-feature policy.")
  }
  adjudication <- NULL
  adjudication_file <- file.path(dirname(mapping_file), "link_adjudication.csv")
  resolved <- retained & nzchar(mapping$conflicted_features)
  if (any(resolved)) {
    if (!all(c("identity_decision", "rejected_features") %in% names(mapping)) ||
        !file.exists(adjudication_file)) stop("Resolved alternatives require explicit link-level adjudication.")
    adjudication <- read.csv(adjudication_file, stringsAsFactors = FALSE)
    for (i in which(resolved)) {
      z <- mapping[i, ]
      links <- adjudication[adjudication$human_symbol == z$human_symbol & adjudication$target_id == z$target_ids, , drop = FALSE]
      rejected <- strsplit(z$conflicted_features, ";", fixed = TRUE)[[1]]
      if (z$identity_decision != "selected_stable_target_resolves_rejected_alternative" ||
          !identical(z$rejected_features, z$conflicted_features) ||
          !identical(links$feature[links$decision == "accepted"], z$feature) ||
          !setequal(links$feature[links$decision == "rejected_alternative"], rejected)) {
        stop("Retained identity does not match explicit link adjudication: ", z$human_symbol)
      }
    }
  }
  if (any(retained & (!nzchar(mapping$feature) | mapping$reason != "retained"))) stop("Invalid approved mapping.")
  mapping$in_universe <- retained & mapping$feature %in% features
  mapping$effective_reason <- mapping$reason
  mapping$effective_reason[retained & !mapping$in_universe] <- "not_in_common_count_universe"
  mapping$score_name <- paste(mapping$signature, mapping$arm, sep = "__")
  arms <- unique(mapping$score_name)
  gene_sets <- setNames(lapply(arms, function(a) sort(unique(mapping$feature[
    mapping$score_name == a & mapping$in_universe]), method = "radix")), arms)
  coverage <- do.call(rbind, lapply(arms, function(a) {
    z <- mapping[mapping$score_name == a, , drop = FALSE]
    n <- length(gene_sets[[a]])
    data.frame(score_name = a, signature = z$signature[1], arm = z$arm[1],
      human_total = nrow(z), human_mapping_retained = sum(z$approved),
      human_in_universe = sum(z$in_universe), human_coverage = mean(z$in_universe),
      effective_pig_genes = n, status = if (n) "scored" else "empty_not_scored",
      small_set_warning = n > 0 && n < 10, stringsAsFactors = FALSE)
  }))
  list(mapping = mapping, gene_sets = gene_sets, coverage = coverage,
    adjudication = adjudication, adjudication_file = adjudication_file)
}
