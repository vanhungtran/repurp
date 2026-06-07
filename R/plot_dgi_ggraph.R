#' Plot a tripartite Drug–Gene–Pathway network with ggraph (full view)
#'
#' Builds a Fruchterman-Reingold `ggraph` plot with three node types:
#' Gene (green), Drug (yellow), and optionally Pathway (purple).
#' Hub drugs (those targeting ≥ `min_drug_targets` genes) are labelled by default.
#'
#' @param dgi_edges A data frame with columns `from` (gene), `to` (drug), `edge_type`
#'   (interaction type string). Typically obtained by calling
#'   `dplyr::select(net_dat, from = Gene, to = Drug, edge_type = InteractionType)`.
#' @param pathway_edges Optional data frame with columns `from` (gene), `to` (pathway
#'   description), `edge_type` (pathway category, e.g. "KEGG"). Pass `NULL` (default)
#'   to omit pathways.
#' @param title Plot title. Default `"Drug–Gene Network"`.
#' @param subtitle Optional subtitle string. Auto-generated when `NULL` (default).
#' @param caption Optional caption. Auto-generated when `NULL`.
#' @param min_drug_targets Minimum number of gene targets for a drug to be labelled
#'   in the full view. Default `2`.
#' @param seed Random seed for reproducible layout. Default `42`.
#' @param layout ggraph layout algorithm. Default `"fr"` (Fruchterman-Reingold).
#' @param out_dir Directory for saved files. `NULL` = plot to current device only.
#' @param base_filename Base filename without extension. Default `paste0("DGI_ggraph_", Sys.Date())`.
#' @param width,height Plot dimensions in inches. Defaults `20` / `18`.
#' @param dpi Resolution for PNG. Default `300`.
#'
#' @return Invisibly returns the `ggplot` object.
#' @export
#' @examples
#' \dontrun{
#' dgi <- dplyr::select(net_dat, from = Gene, to = Drug, edge_type = InteractionType)
#' plot_dgi_ggraph(dgi, title = "My network", out_dir = "results/")
#' }
plot_dgi_ggraph <- function(dgi_edges,
                             pathway_edges    = NULL,
                             title            = "Drug–Gene Network",
                             subtitle         = NULL,
                             caption          = NULL,
                             min_drug_targets = 2L,
                             seed             = 42L,
                             layout           = "fr",
                             out_dir          = NULL,
                             base_filename    = paste0("DGI_ggraph_", Sys.Date()),
                             width  = 20,
                             height = 18,
                             dpi    = 300) {

  .check_ggraph()

  dgi_edges$edge_type <- .normalise_interaction(dgi_edges$edge_type)

  gene_ids <- unique(dgi_edges$from)
  drug_ids <- unique(dgi_edges$to)
  drug_degree <- dplyr::count(dgi_edges, to, name = "targets")

  nodes <- dplyr::bind_rows(
    tibble::tibble(name = gene_ids, type = "Gene"),
    tibble::tibble(name = drug_ids, type = "Drug") |>
      dplyr::left_join(drug_degree, by = c("name" = "to")) |>
      dplyr::mutate(targets = ifelse(is.na(targets), 0L, targets))
  )

  all_edges <- dgi_edges

  if (!is.null(pathway_edges)) {
    path_ids <- unique(pathway_edges$to)
    nodes <- dplyr::bind_rows(
      nodes,
      tibble::tibble(name = path_ids, type = "Pathway", targets = NA_integer_)
    )
    all_edges <- dplyr::bind_rows(dgi_edges, pathway_edges)
  }

  nodes <- dplyr::mutate(nodes, id = dplyr::row_number())

  edge_df <- all_edges |>
    dplyr::left_join(nodes, by = c("from" = "name")) |>
    dplyr::rename(from_id = id) |>
    dplyr::left_join(nodes, by = c("to" = "name")) |>
    dplyr::rename(to_id = id) |>
    dplyr::select(from = from_id, to = to_id, edge_type)

  graph <- tidygraph::tbl_graph(nodes = nodes, edges = edge_df, directed = FALSE)

  node_pal <- dgi_node_colors
  edge_pal <- c(dgi_edge_colors,
                KEGG = "#1F78B4", WikiPathways = "#33A02C", GO_Process = "#FB9A99")

  graph <- graph |>
    tidygraph::activate("nodes") |>
    dplyr::mutate(
      show_label = dplyr::case_when(
        type == "Gene"    ~ TRUE,
        type == "Pathway" ~ TRUE,
        type == "Drug"    ~ (!is.na(targets) & targets >= min_drug_targets)
      ),
      label_size = dplyr::case_when(
        type == "Gene"    ~ 2.5,
        type == "Pathway" ~ 3.0,
        type == "Drug"    ~ 2.0
      ),
      label_face = dplyr::case_when(
        type == "Gene"    ~ "plain",
        type == "Pathway" ~ "bold",
        type == "Drug"    ~ "italic"
      ),
      point_size = dplyr::case_when(
        type == "Gene"    ~ 3.0,
        type == "Drug"    ~ ifelse(!is.na(targets) & targets >= min_drug_targets,
                                   2.5, 0.8),
        type == "Pathway" ~ 5.0
      )
    )

  if (is.null(subtitle)) {
    n_hub <- sum(!is.na(drug_degree$targets) & drug_degree$targets >= min_drug_targets)
    subtitle <- paste0(
      length(gene_ids), " genes | ",
      length(drug_ids), " drugs (", n_hub, " hub drugs labelled)",
      if (!is.null(pathway_edges)) paste0(" | ", length(unique(pathway_edges$to)), " pathways") else ""
    )
  }
  if (is.null(caption)) {
    caption <- paste0("DGIdb × Repurp | ", Sys.Date())
  }

  set.seed(seed)
  p <- ggraph::ggraph(graph, layout = layout) +
    ggraph::geom_edge_link(ggplot2::aes(color = edge_type),
                            alpha = 0.2, width = 0.25) +
    ggraph::geom_node_point(ggplot2::aes(color = type, size = point_size),
                             alpha = 0.85) +
    ggraph::geom_node_text(
      data = function(d) d[d$show_label, ],
      ggplot2::aes(label = name, size = label_size, fontface = label_face),
      repel = TRUE, max.overlaps = 60, family = "sans",
      segment.size = 0.15, segment.alpha = 0.4,
      box.padding = 0.3, force = 0.5
    ) +
    ggraph::scale_edge_colour_manual(values = edge_pal, name = NULL) +
    ggplot2::scale_color_manual(values = node_pal, name = NULL) +
    ggplot2::scale_size_identity() +
    ggplot2::guides(color = "none", edge_colour = "none") +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption) +
    ggraph::theme_graph(background = "white") +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5, color = "grey40"),
      plot.caption  = ggplot2::element_text(size = 7, color = "grey60"),
      plot.margin   = ggplot2::margin(15, 15, 15, 15)
    )

  .save_ggplot(p, out_dir, base_filename, width, height, dpi)
  invisible(p)
}


#' Plot a core Drug–Gene–Pathway network (high-connectivity nodes only)
#'
#' Filters to genes with at least `min_gene_edges` edges before plotting,
#' which reduces clutter for large networks. All drug and pathway labels are shown.
#'
#' @inheritParams plot_dgi_ggraph
#' @param min_gene_edges Minimum edge count for a gene to be included. Default `3`.
#'
#' @return Invisibly returns the `ggplot` object.
#' @export
#' @examples
#' \dontrun{
#' dgi <- dplyr::select(net_dat, from = Gene, to = Drug, edge_type = InteractionType)
#' plot_dgi_ggraph_core(dgi, out_dir = "results/", min_gene_edges = 3)
#' }
plot_dgi_ggraph_core <- function(dgi_edges,
                                  pathway_edges   = NULL,
                                  title           = "Core Drug–Gene Network",
                                  subtitle        = "Genes with ≥ 3 edges | All drugs and pathways labelled",
                                  caption         = NULL,
                                  min_gene_edges  = 3L,
                                  seed            = 42L,
                                  layout          = "kk",
                                  out_dir         = NULL,
                                  base_filename   = paste0("DGI_ggraph_core_", Sys.Date()),
                                  width  = 16,
                                  height = 14,
                                  dpi    = 300) {

  .check_ggraph()

  dgi_edges$edge_type <- .normalise_interaction(dgi_edges$edge_type)

  gene_ids <- unique(dgi_edges$from)
  drug_ids <- unique(dgi_edges$to)
  drug_degree <- dplyr::count(dgi_edges, to, name = "targets")

  nodes <- dplyr::bind_rows(
    tibble::tibble(name = gene_ids, type = "Gene"),
    tibble::tibble(name = drug_ids, type = "Drug") |>
      dplyr::left_join(drug_degree, by = c("name" = "to")) |>
      dplyr::mutate(targets = ifelse(is.na(targets), 0L, targets))
  )

  all_edges <- dgi_edges
  if (!is.null(pathway_edges)) {
    path_ids <- unique(pathway_edges$to)
    nodes <- dplyr::bind_rows(
      nodes,
      tibble::tibble(name = path_ids, type = "Pathway", targets = NA_integer_)
    )
    all_edges <- dplyr::bind_rows(dgi_edges, pathway_edges)
  }

  nodes <- dplyr::mutate(nodes, id = dplyr::row_number())

  edge_df <- all_edges |>
    dplyr::left_join(nodes, by = c("from" = "name")) |>
    dplyr::rename(from_id = id) |>
    dplyr::left_join(nodes, by = c("to" = "name")) |>
    dplyr::rename(to_id = id) |>
    dplyr::select(from = from_id, to = to_id, edge_type)

  graph <- tidygraph::tbl_graph(nodes = nodes, edges = edge_df, directed = FALSE)

  # Filter to core genes
  gene_edge_counts <- dplyr::count(edge_df, from, name = "n")
  core_gene_ids <- nodes$id[nodes$type == "Gene" &
                              nodes$id %in% gene_edge_counts$from[
                                gene_edge_counts$n >= min_gene_edges]]

  graph <- graph |>
    tidygraph::activate("nodes") |>
    dplyr::filter(id %in% core_gene_ids | type != "Gene") |>
    dplyr::mutate(deg = tidygraph::centrality_degree()) |>
    dplyr::filter(deg > 0) |>
    dplyr::mutate(
      show_label = TRUE,
      label_size = dplyr::case_when(
        type == "Gene"    ~ 3.0,
        type == "Pathway" ~ 3.5,
        type == "Drug"    ~ 2.2
      ),
      label_face = dplyr::case_when(
        type == "Gene"    ~ "bold",
        type == "Pathway" ~ "bold",
        type == "Drug"    ~ "plain"
      ),
      point_size = dplyr::case_when(
        type == "Gene"    ~ 4.0,
        type == "Drug"    ~ 2.5,
        type == "Pathway" ~ 7.0
      )
    )

  if (is.null(caption)) caption <- paste0("DGIdb × Repurp | ", Sys.Date())

  node_pal <- dgi_node_colors
  edge_pal <- c(dgi_edge_colors,
                KEGG = "#1F78B4", WikiPathways = "#33A02C", GO_Process = "#FB9A99")

  set.seed(seed)
  p <- ggraph::ggraph(graph, layout = layout) +
    ggraph::geom_edge_link(ggplot2::aes(color = edge_type),
                            alpha = 0.3, width = 0.5) +
    ggraph::geom_node_point(ggplot2::aes(color = type, size = point_size),
                             alpha = 0.9) +
    ggraph::geom_node_text(
      ggplot2::aes(label = name, size = label_size, fontface = label_face),
      repel = TRUE, max.overlaps = 80, family = "sans",
      segment.size = 0.2, segment.alpha = 0.5,
      box.padding = 0.4, force = 2
    ) +
    ggraph::scale_edge_colour_manual(values = edge_pal, name = NULL) +
    ggplot2::scale_color_manual(values = node_pal, name = NULL) +
    ggplot2::scale_size_identity() +
    ggplot2::guides(color = "none", edge_colour = "none") +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption) +
    ggraph::theme_graph(background = "white") +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5, color = "grey40"),
      plot.caption  = ggplot2::element_text(size = 7, color = "grey60"),
      plot.margin   = ggplot2::margin(15, 15, 15, 15)
    )

  .save_ggplot(p, out_dir, base_filename, width, height, dpi)
  invisible(p)
}
