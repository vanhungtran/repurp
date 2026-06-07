#' Prepare a drug-gene interaction data frame for network plotting
#'
#' Reads a DGIdb-style interactions TSV (or accepts a pre-loaded data frame) and
#' filters to a user-supplied gene list, returning a standardised four-column
#' data frame ready for all `plot_dgi_*` functions.
#'
#' @param interactions A data frame with at least columns `gene_name`, `drug_name`,
#'   and `interaction_type` (DGIdb TSV format), **or** the path to such a TSV file.
#' @param gene_list Character vector of gene symbols to retain.
#' @param source_col Optional name of a source-database column to keep (string).
#'   Defaults to `"interaction_source_db_name"` when present, otherwise omitted.
#' @param interaction_score_col Optional column name for a numeric interaction score.
#'   Kept as-is when present.
#'
#' @return A tibble with columns `Gene`, `Drug`, `InteractionType`, and optionally
#'   `Source` and `interaction_score`.
#'
#' @export
#' @examples
#' \dontrun{
#' net_dat <- prepare_dgi_data("interactions.tsv", gene_list = my_genes)
#' }
prepare_dgi_data <- function(interactions, gene_list,
                              source_col = NULL,
                              interaction_score_col = NULL) {
  if (is.character(interactions) && length(interactions) == 1 &&
      file.exists(interactions)) {
    interactions <- utils::read.delim(interactions, stringsAsFactors = FALSE)
  }

  if (!inherits(interactions, "data.frame")) {
    rlang::abort("`interactions` must be a data frame or a valid file path.")
  }

  required <- c("gene_name", "drug_name", "interaction_type")
  missing_cols <- setdiff(required, colnames(interactions))
  if (length(missing_cols) > 0) {
    rlang::abort(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  # Detect source column automatically
  if (is.null(source_col) && "interaction_source_db_name" %in% colnames(interactions)) {
    source_col <- "interaction_source_db_name"
  }

  keep_cols <- c("gene_name", "drug_name", "interaction_type")
  if (!is.null(source_col) && source_col %in% colnames(interactions)) {
    keep_cols <- c(keep_cols, source_col)
  }
  if (!is.null(interaction_score_col) &&
      interaction_score_col %in% colnames(interactions)) {
    keep_cols <- c(keep_cols, interaction_score_col)
  }

  out <- interactions[interactions$gene_name %in% gene_list, keep_cols,
                      drop = FALSE]

  # Rename standard columns
  colnames(out)[colnames(out) == "gene_name"]        <- "Gene"
  colnames(out)[colnames(out) == "drug_name"]        <- "Drug"
  colnames(out)[colnames(out) == "interaction_type"] <- "InteractionType"
  if (!is.null(source_col) && source_col %in% colnames(out)) {
    colnames(out)[colnames(out) == source_col] <- "Source"
  }

  out$InteractionType <- .normalise_interaction(out$InteractionType)

  # Disambiguate drug names that collide with gene names
  collide <- out$Drug %in% out$Gene
  out$Drug[collide] <- paste0(out$Drug[collide], "_drug")

  tibble::as_tibble(out)
}
