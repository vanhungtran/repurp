#' Standardise drug-gene interaction edge columns
#'
#' Converts common DGIdb and Repurp edge-list layouts into a consistent tibble
#' with `gene`, `drug`, `interaction_type`, and optional `source` columns. This
#' is useful before ranking candidate drugs or combining interaction tables from
#' different sources.
#'
#' @param edges A data frame containing drug-gene interactions.
#' @param gene_col,drug_col,interaction_col,source_col Optional explicit column
#'   names. When omitted, common column names are detected automatically.
#' @param direction Direction of `from`/`to` edge lists. `"gene_to_drug"` treats
#'   `from` as gene and `to` as drug; `"drug_to_gene"` reverses them; `"auto"`
#'   uses `drug_info` overlap when possible and otherwise assumes
#'   gene-to-drug.
#' @param drug_info Optional drug metadata with a `drug_name` column, used only
#'   for automatic direction detection.
#'
#' @return A tibble with columns `gene`, `drug`, `interaction_type`, and
#'   optional `source`.
#' @export
#' @examples
#' data(example_dgi)
#' standardize_dgi_edges(example_dgi)
standardize_dgi_edges <- function(edges,
                                  gene_col = NULL,
                                  drug_col = NULL,
                                  interaction_col = NULL,
                                  source_col = NULL,
                                  direction = c("auto", "gene_to_drug", "drug_to_gene"),
                                  drug_info = repurp_ad_drugs()) {
  direction <- match.arg(direction)
  .standardize_dgi_edges(
    edges = edges,
    gene_col = gene_col,
    drug_col = drug_col,
    interaction_col = interaction_col,
    source_col = source_col,
    direction = direction,
    drug_info = drug_info
  )
}

#' Rank candidate drugs by biomarker coverage and evidence weight
#'
#' Scores drugs by the weighted biomarker targets they cover. Use this after
#' [prepare_dgi_data()] or with any edge list accepted by
#' [standardize_dgi_edges()].
#'
#' @param edges Drug-gene interaction data frame.
#' @param biomarkers Optional character vector defining the biomarker universe
#'   to score against. When omitted, all genes in `edges` are used.
#' @param gene_weights Optional named numeric vector, or a data frame with gene
#'   and weight columns, giving biomarker importance weights. Missing genes use
#'   weight `1`.
#' @param drug_info Optional data frame with `drug_name`, `category`, `brand`,
#'   and/or `targets_known` columns to attach to the output.
#' @param direction Direction for `from`/`to` edge lists; see
#'   [standardize_dgi_edges()].
#' @param interaction_weights Optional named numeric vector weighting
#'   interaction types.
#' @param source_weights Optional named numeric vector weighting evidence
#'   sources.
#' @param direct_bonus Multiplicative bonus for edges whose interaction type is
#'   `"direct"`. Defaults to `1.5`.
#' @param min_targets Minimum number of unique targeted biomarkers required in
#'   the returned table.
#'
#' @return A tibble ranked by descending `score`.
#' @export
#' @examples
#' data(example_dgi)
#' score_drug_repurposing(example_dgi, biomarkers = c("IL13", "IL4R", "JAK1"))
score_drug_repurposing <- function(edges,
                                   biomarkers = NULL,
                                   gene_weights = NULL,
                                   drug_info = repurp_ad_drugs(),
                                   direction = c("auto", "gene_to_drug", "drug_to_gene"),
                                   interaction_weights = NULL,
                                   source_weights = NULL,
                                   direct_bonus = 1.5,
                                   min_targets = 1L) {
  direction <- match.arg(direction)
  dgi <- .standardize_dgi_edges(edges, direction = direction, drug_info = drug_info)

  if (!is.null(biomarkers)) {
    biomarkers <- unique(as.character(biomarkers))
    dgi <- dgi[dgi$gene %in% biomarkers, , drop = FALSE]
    universe <- biomarkers
  } else {
    universe <- unique(dgi$gene)
  }

  if (nrow(dgi) == 0) {
    return(.empty_repurposing_score())
  }

  gene_weight_tbl <- .coerce_gene_weights(gene_weights, universe)
  dgi <- dplyr::left_join(dgi, gene_weight_tbl, by = "gene")
  dgi$gene_weight[is.na(dgi$gene_weight)] <- 1

  dgi$interaction_weight <- .lookup_named_weight(
    dgi$interaction_type,
    interaction_weights,
    default = 1
  )
  dgi$source_weight <- if ("source" %in% colnames(dgi)) {
    .lookup_named_weight(dgi$source, source_weights, default = 1)
  } else {
    1
  }
  dgi$direct_multiplier <- ifelse(
    tolower(dgi$interaction_type) == "direct",
    direct_bonus,
    1
  )
  dgi$edge_score <- dgi$gene_weight * dgi$interaction_weight *
    dgi$source_weight * dgi$direct_multiplier

  per_gene <- dgi |>
    dplyr::group_by(drug, gene) |>
    dplyr::summarise(
      gene_score = max(edge_score, na.rm = TRUE),
      interactions = paste(sort(unique(interaction_type)), collapse = "; "),
      sources = if ("source" %in% colnames(dgi)) {
        .collapse_terms(source)
      } else {
        ""
      },
      .groups = "drop"
    )

  out <- per_gene |>
    dplyr::group_by(drug) |>
    dplyr::summarise(
      score = sum(gene_score, na.rm = TRUE),
      target_count = dplyr::n_distinct(gene),
      coverage = target_count / length(universe),
      targets = paste(sort(unique(gene)), collapse = "; "),
      interaction_types = .collapse_terms(interactions),
      sources = .collapse_terms(sources),
      .groups = "drop"
    ) |>
    dplyr::filter(target_count >= min_targets) |>
    dplyr::arrange(dplyr::desc(score), dplyr::desc(target_count), drug)

  out$rank <- seq_len(nrow(out))
  out <- out[, c("rank", setdiff(colnames(out), "rank")), drop = FALSE]

  .attach_drug_info(out, drug_info)
}

#' Summarise biomarker coverage by available drugs
#'
#' Reports how many drugs and drug classes cover each biomarker, including
#' uncovered biomarkers when `biomarkers` is supplied.
#'
#' @param edges Drug-gene interaction data frame.
#' @param biomarkers Optional biomarker universe to include in the output.
#' @param drug_info Optional data frame with `drug_name` and `category` columns.
#' @param direction Direction for `from`/`to` edge lists; see
#'   [standardize_dgi_edges()].
#'
#' @return A tibble with one row per biomarker.
#' @export
#' @examples
#' data(example_dgi)
#' summarize_biomarker_coverage(example_dgi)
summarize_biomarker_coverage <- function(edges,
                                         biomarkers = NULL,
                                         drug_info = repurp_ad_drugs(),
                                         direction = c("auto", "gene_to_drug", "drug_to_gene")) {
  direction <- match.arg(direction)
  dgi <- .standardize_dgi_edges(edges, direction = direction, drug_info = drug_info)

  if (!is.null(biomarkers)) {
    biomarkers <- unique(as.character(biomarkers))
    dgi <- dgi[dgi$gene %in% biomarkers, , drop = FALSE]
    universe <- tibble::tibble(gene = biomarkers)
  } else {
    universe <- tibble::tibble(gene = sort(unique(dgi$gene)))
  }

  dgi <- .add_drug_category(dgi, drug_info)

  summary <- dgi |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      drug_count = dplyr::n_distinct(drug),
      drug_class_count = dplyr::n_distinct(category[!is.na(category)]),
      direct_drug_count = dplyr::n_distinct(drug[tolower(interaction_type) == "direct"]),
      drugs = paste(sort(unique(drug)), collapse = "; "),
      drug_classes = paste(sort(unique(category[!is.na(category)])), collapse = "; "),
      interaction_types = paste(sort(unique(interaction_type)), collapse = "; "),
      .groups = "drop"
    )

  out <- dplyr::left_join(universe, summary, by = "gene")
  count_cols <- c("drug_count", "drug_class_count", "direct_drug_count")
  for (col in count_cols) {
    out[[col]][is.na(out[[col]])] <- 0L
  }
  text_cols <- c("drugs", "drug_classes", "interaction_types")
  for (col in text_cols) {
    out[[col]][is.na(out[[col]])] <- ""
  }

  out |>
    dplyr::arrange(dplyr::desc(drug_count), gene)
}

#' Summarise drug classes represented in a repurposing network
#'
#' Aggregates drug-gene edges by pharmacological class, reporting how many drugs
#' and biomarkers each class covers.
#'
#' @param edges Drug-gene interaction data frame.
#' @param biomarkers Optional biomarker universe used to calculate coverage.
#' @param drug_info Drug metadata with `drug_name` and `category` columns.
#' @param direction Direction for `from`/`to` edge lists; see
#'   [standardize_dgi_edges()].
#'
#' @return A tibble ranked by descending biomarker coverage and edge count.
#' @export
#' @examples
#' ad_edges <- repurp_ad_edges()
#' summarize_drug_classes(ad_edges, direction = "drug_to_gene")
summarize_drug_classes <- function(edges,
                                   biomarkers = NULL,
                                   drug_info = repurp_ad_drugs(),
                                   direction = c("auto", "gene_to_drug", "drug_to_gene")) {
  direction <- match.arg(direction)
  dgi <- .standardize_dgi_edges(edges, direction = direction, drug_info = drug_info)

  if (!is.null(biomarkers)) {
    biomarkers <- unique(as.character(biomarkers))
    dgi <- dgi[dgi$gene %in% biomarkers, , drop = FALSE]
    denominator <- length(biomarkers)
  } else {
    denominator <- dplyr::n_distinct(dgi$gene)
  }

  if (nrow(dgi) == 0) {
    return(tibble::tibble(
      category = character(),
      drug_count = integer(),
      target_count = integer(),
      coverage = numeric(),
      edge_count = integer(),
      drugs = character(),
      targets = character()
    ))
  }

  dgi <- .add_drug_category(dgi, drug_info)
  dgi$category[is.na(dgi$category) | dgi$category == ""] <- "Other"

  dgi |>
    dplyr::group_by(category) |>
    dplyr::summarise(
      drug_count = dplyr::n_distinct(drug),
      target_count = dplyr::n_distinct(gene),
      coverage = target_count / denominator,
      edge_count = dplyr::n(),
      drugs = paste(sort(unique(drug)), collapse = "; "),
      targets = paste(sort(unique(gene)), collapse = "; "),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(coverage), dplyr::desc(edge_count), category)
}

#' Score drugs by expected reversal of disease gene expression
#'
#' Combines drug-gene interactions with differential-expression results. A drug
#' receives a positive contribution when its interaction type is expected to
#' oppose the disease direction of a target gene, for example inhibition of an
#' up-regulated biomarker or activation of a down-regulated biomarker.
#'
#' @param edges Drug-gene interaction data frame.
#' @param gene_stats Data frame containing gene symbols and log fold-changes.
#' @param gene_col,logfc_col,p_col Optional explicit column names in
#'   `gene_stats`. `p_col` is optional; when present, `-log10(p)` weights the
#'   expression signal.
#' @param drug_info Optional drug metadata attached to the output.
#' @param direction Direction for `from`/`to` edge lists; see
#'   [standardize_dgi_edges()].
#' @param interaction_effects Optional named numeric vector mapping interaction
#'   types to expected drug effect on the target: `-1` decreases/inhibits,
#'   `1` increases/activates, `0` unknown.
#' @param min_directional_edges Minimum number of edges with known effect
#'   direction required in the returned table.
#'
#' @return A tibble ranked by descending `reversal_score`.
#' @export
#' @examples
#' data(example_dgi)
#' deg <- tibble::tibble(gene = c("IL13", "IL4R"), logFC = c(1.2, -0.8))
#' score_expression_reversal(example_dgi, deg)
score_expression_reversal <- function(edges,
                                      gene_stats,
                                      gene_col = NULL,
                                      logfc_col = NULL,
                                      p_col = NULL,
                                      drug_info = repurp_ad_drugs(),
                                      direction = c("auto", "gene_to_drug", "drug_to_gene"),
                                      interaction_effects = NULL,
                                      min_directional_edges = 1L) {
  direction <- match.arg(direction)
  dgi <- .standardize_dgi_edges(edges, direction = direction, drug_info = drug_info)
  stats <- .standardize_gene_stats(gene_stats, gene_col, logfc_col, p_col)

  dgi <- dplyr::inner_join(dgi, stats, by = "gene")
  if (nrow(dgi) == 0) {
    return(.empty_reversal_score())
  }

  dgi$effect_direction <- .interaction_effect_direction(
    dgi$interaction_type,
    interaction_effects
  )
  dgi$expression_weight <- abs(dgi$logFC)
  if ("p_value" %in% colnames(dgi)) {
    clean_p <- pmax(as.numeric(dgi$p_value), .Machine$double.xmin, na.rm = TRUE)
    dgi$expression_weight <- dgi$expression_weight * -log10(clean_p)
  }
  dgi$edge_reversal_score <- -sign(dgi$logFC) * dgi$effect_direction *
    dgi$expression_weight

  out <- dgi |>
    dplyr::group_by(drug) |>
    dplyr::summarise(
      reversal_score = sum(edge_reversal_score, na.rm = TRUE),
      mean_edge_score = mean(edge_reversal_score, na.rm = TRUE),
      directional_edges = sum(effect_direction != 0, na.rm = TRUE),
      supportive_edges = sum(edge_reversal_score > 0, na.rm = TRUE),
      conflicting_edges = sum(edge_reversal_score < 0, na.rm = TRUE),
      unknown_edges = sum(effect_direction == 0 | is.na(effect_direction), na.rm = TRUE),
      target_count = dplyr::n_distinct(gene),
      targets = paste(sort(unique(gene)), collapse = "; "),
      interaction_types = paste(sort(unique(interaction_type)), collapse = "; "),
      .groups = "drop"
    ) |>
    dplyr::filter(directional_edges >= min_directional_edges) |>
    dplyr::arrange(dplyr::desc(reversal_score), dplyr::desc(supportive_edges), drug)

  out$rank <- seq_len(nrow(out))
  out <- out[, c("rank", setdiff(colnames(out), "rank")), drop = FALSE]

  .attach_drug_info(out, drug_info)
}

.standardize_dgi_edges <- function(edges,
                                   gene_col = NULL,
                                   drug_col = NULL,
                                   interaction_col = NULL,
                                   source_col = NULL,
                                   direction = c("auto", "gene_to_drug", "drug_to_gene"),
                                   drug_info = repurp_ad_drugs()) {
  direction <- match.arg(direction)

  if (!inherits(edges, "data.frame")) {
    rlang::abort("`edges` must be a data frame.")
  }
  if (nrow(edges) == 0) {
    rlang::abort("`edges` has no rows.")
  }

  nms <- colnames(edges)
  has_from_to <- all(c("from", "to") %in% nms)

  if (is.null(gene_col) && is.null(drug_col) && has_from_to &&
      !any(c("Gene", "gene", "gene_name") %in% nms) &&
      !any(c("Drug", "drug", "drug_name") %in% nms)) {
    resolved <- .resolve_edge_direction(edges, direction, drug_info)
    gene_col <- if (resolved == "drug_to_gene") "to" else "from"
    drug_col <- if (resolved == "drug_to_gene") "from" else "to"
  } else {
    gene_col <- .detect_column(
      edges,
      explicit = gene_col,
      candidates = c("gene", "Gene", "gene_name", "symbol", "Symbol"),
      label = "gene"
    )
    drug_col <- .detect_column(
      edges,
      explicit = drug_col,
      candidates = c("drug", "Drug", "drug_name", "compound", "Compound"),
      label = "drug"
    )
  }

  interaction_col <- .detect_column(
    edges,
    explicit = interaction_col,
    candidates = c("interaction_type", "InteractionType", "edge_type", "type"),
    label = "interaction",
    required = FALSE
  )
  source_col <- .detect_column(
    edges,
    explicit = source_col,
    candidates = c("source", "Source", "interaction_source_db_name", "database"),
    label = "source",
    required = FALSE
  )

  out <- tibble::tibble(
    gene = as.character(edges[[gene_col]]),
    drug = as.character(edges[[drug_col]]),
    interaction_type = if (!is.null(interaction_col)) {
      as.character(edges[[interaction_col]])
    } else {
      "Other"
    }
  )
  out$interaction_type <- .normalise_interaction(out$interaction_type)

  if (!is.null(source_col)) {
    out$source <- as.character(edges[[source_col]])
  }

  out <- out[!is.na(out$gene) & out$gene != "" &
               !is.na(out$drug) & out$drug != "", , drop = FALSE]
  dplyr::distinct(out)
}

.detect_column <- function(data, explicit, candidates, label, required = TRUE) {
  if (!is.null(explicit)) {
    if (!explicit %in% colnames(data)) {
      rlang::abort(sprintf("Column `%s` was not found for %s.", explicit, label))
    }
    return(explicit)
  }

  hit <- candidates[candidates %in% colnames(data)]
  if (length(hit) > 0) {
    return(hit[[1]])
  }

  if (required) {
    rlang::abort(sprintf(
      "Could not detect a %s column. Supply `%s_col` explicitly.",
      label,
      label
    ))
  }
  NULL
}

.resolve_edge_direction <- function(edges, direction, drug_info) {
  if (direction != "auto") {
    return(direction)
  }

  known_drugs <- character()
  if (inherits(drug_info, "data.frame") && "drug_name" %in% colnames(drug_info)) {
    known_drugs <- as.character(drug_info$drug_name)
  }

  from_hits <- sum(as.character(edges$from) %in% known_drugs, na.rm = TRUE)
  to_hits <- sum(as.character(edges$to) %in% known_drugs, na.rm = TRUE)

  if (from_hits > to_hits) {
    "drug_to_gene"
  } else {
    "gene_to_drug"
  }
}

.coerce_gene_weights <- function(gene_weights, universe) {
  if (is.null(gene_weights)) {
    return(tibble::tibble(gene = universe, gene_weight = rep(1, length(universe))))
  }

  if (is.numeric(gene_weights) && !is.null(names(gene_weights))) {
    return(tibble::tibble(
      gene = names(gene_weights),
      gene_weight = as.numeric(gene_weights)
    ))
  }

  if (inherits(gene_weights, "data.frame")) {
    gene_col <- .detect_column(
      gene_weights,
      explicit = NULL,
      candidates = c("gene", "Gene", "gene_name", "symbol", "Symbol"),
      label = "gene"
    )
    weight_col <- .detect_column(
      gene_weights,
      explicit = NULL,
      candidates = c("weight", "Weight", "gene_weight", "importance", "score"),
      label = "weight"
    )
    return(tibble::tibble(
      gene = as.character(gene_weights[[gene_col]]),
      gene_weight = as.numeric(gene_weights[[weight_col]])
    ))
  }

  rlang::abort(
    "`gene_weights` must be NULL, a named numeric vector, or a data frame with gene and weight columns."
  )
}

.lookup_named_weight <- function(values, weights, default = 1) {
  if (is.null(weights)) {
    return(rep(default, length(values)))
  }
  if (!is.numeric(weights) || is.null(names(weights))) {
    rlang::abort("Weight mappings must be named numeric vectors.")
  }

  out <- rep(default, length(values))
  idx <- match(values, names(weights))
  matched <- !is.na(idx)
  out[matched] <- as.numeric(weights[idx[matched]])
  out
}

.collapse_terms <- function(x) {
  if (length(x) == 0) {
    return("")
  }
  parts <- unlist(strsplit(as.character(x), "; ", fixed = TRUE), use.names = FALSE)
  parts <- sort(unique(parts[!is.na(parts) & parts != ""]))
  paste(parts, collapse = "; ")
}

.attach_drug_info <- function(out, drug_info) {
  if (!inherits(drug_info, "data.frame") || !"drug_name" %in% colnames(drug_info)) {
    return(out)
  }

  keep <- intersect(
    c("drug_name", "brand", "category", "targets_known"),
    colnames(drug_info)
  )
  info <- dplyr::distinct(drug_info[, keep, drop = FALSE])
  out <- dplyr::left_join(out, info, by = c("drug" = "drug_name"))

  front <- intersect(c("rank", "drug", "brand", "category", "targets_known"), colnames(out))
  out[, c(front, setdiff(colnames(out), front)), drop = FALSE]
}

.add_drug_category <- function(dgi, drug_info) {
  if (inherits(drug_info, "data.frame") &&
      all(c("drug_name", "category") %in% colnames(drug_info))) {
    dgi <- dplyr::left_join(
      dgi,
      dplyr::distinct(drug_info[, c("drug_name", "category"), drop = FALSE]),
      by = c("drug" = "drug_name")
    )
  } else {
    dgi$category <- NA_character_
  }
  dgi
}

.standardize_gene_stats <- function(gene_stats, gene_col, logfc_col, p_col) {
  if (!inherits(gene_stats, "data.frame")) {
    rlang::abort("`gene_stats` must be a data frame.")
  }

  gene_col <- .detect_column(
    gene_stats,
    explicit = gene_col,
    candidates = c("gene", "Gene", "gene_name", "symbol", "Symbol"),
    label = "gene"
  )
  logfc_col <- .detect_column(
    gene_stats,
    explicit = logfc_col,
    candidates = c("logFC", "log2FoldChange", "avg_log2FC", "estimate", "effect"),
    label = "logfc"
  )
  p_col <- .detect_column(
    gene_stats,
    explicit = p_col,
    candidates = c("adj.P.Val", "padj", "FDR", "p_val_adj", "pvalue", "p.value", "P.Value"),
    label = "p",
    required = FALSE
  )

  out <- tibble::tibble(
    gene = as.character(gene_stats[[gene_col]]),
    logFC = as.numeric(gene_stats[[logfc_col]])
  )
  if (!is.null(p_col)) {
    out$p_value <- as.numeric(gene_stats[[p_col]])
  }

  out[!is.na(out$gene) & out$gene != "" & !is.na(out$logFC), , drop = FALSE]
}

.interaction_effect_direction <- function(interaction_type, interaction_effects) {
  defaults <- c(
    INHIBITOR = -1,
    ANTAGONIST = -1,
    BLOCKER = -1,
    Antibody = -1,
    Cleavage = -1,
    `Inverse Agonist` = -1,
    AGONIST = 1,
    ACTIVATOR = 1,
    SUBSTRATE = 1,
    MODULATOR = 0,
    BINDER = 0,
    Other = 0,
    None = 0,
    direct = 0,
    indirect = 0
  )

  if (!is.null(interaction_effects)) {
    if (!is.numeric(interaction_effects) || is.null(names(interaction_effects))) {
      rlang::abort("`interaction_effects` must be a named numeric vector.")
    }
    defaults[names(interaction_effects)] <- interaction_effects
  }

  idx <- match(tolower(interaction_type), tolower(names(defaults)))
  out <- rep(0, length(interaction_type))
  matched <- !is.na(idx)
  out[matched] <- as.numeric(defaults[idx[matched]])
  out
}

.empty_repurposing_score <- function() {
  tibble::tibble(
    rank = integer(),
    drug = character(),
    score = numeric(),
    target_count = integer(),
    coverage = numeric(),
    targets = character(),
    interaction_types = character(),
    sources = character()
  )
}

.empty_reversal_score <- function() {
  tibble::tibble(
    rank = integer(),
    drug = character(),
    reversal_score = numeric(),
    mean_edge_score = numeric(),
    directional_edges = integer(),
    supportive_edges = integer(),
    conflicting_edges = integer(),
    unknown_edges = integer(),
    target_count = integer(),
    targets = character(),
    interaction_types = character()
  )
}
