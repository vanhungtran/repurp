#' Plot an enhanced Drug–Gene network with ggraph + qgraph FR layout + ggrepel
#'
#' Builds a publication-quality bipartite drug–gene network using a manual
#' Fruchterman-Reingold layout computed by `qgraph`, rendered with `ggraph`,
#' `ggrepel` labels, and `ggnewscale` for separate node-type and edge-type
#' colour legends.
#'
#' @param dgi_edges A data frame with columns `from` (gene), `to` (drug),
#'   `edge_type` (interaction type string). An optional `score` column
#'   is used for edge width scaling.
#' @param title Plot title.
#' @param mode Either `"full"` (show hub drugs labelled, genes with ≥ `min_gene_edges`)
#'   or `"core"` (filter to core genes, label all nodes). Default `"full"`.
#' @param min_gene_edges Minimum number of drug connections for a gene to be
#'   included in core mode, or labelled in full mode.
#' @param min_drug_targets Minimum number of gene targets for a drug to be labelled
#'   in full mode. Default `2`.
#' @param seed Random seed for reproducible layout.
#' @param out_dir Directory for saved files. `NULL` = plot to current device only.
#' @param base_filename Base filename without extension.
#' @param width,height Plot dimensions in inches. Defaults `20` / `18`.
#' @param dpi Resolution for PNG. Default `300`.
#'
#' @return Invisibly returns the `ggplot` object.
#' @export
#' @examples
#' \dontrun{
#' dgi <- dplyr::transmute(net_dat, from = Gene, to = Drug,
#'                          edge_type = InteractionType, score = interaction_score)
#' plot_dgi_ggraph_enhanced(dgi, mode = "full", out_dir = "results/")
#' plot_dgi_ggraph_enhanced(dgi, mode = "core", out_dir = "results/", min_gene_edges = 3)
#' }
plot_dgi_ggraph_enhanced <- function(dgi_edges,
                                      title            = "Drug–Gene Interaction Network",
                                      mode             = c("full", "core"),
                                      min_gene_edges   = 3L,
                                      min_drug_targets = 2L,
                                      seed             = 42L,
                                      out_dir          = NULL,
                                      base_filename    = paste0("DGI_ggraph_enhanced_", Sys.Date()),
                                      width  = 20,
                                      height = 18,
                                      dpi    = 300) {

  mode <- match.arg(mode)
  .check_ggraph()
  .check_optional_pkg("ggrepel", "ggraph_enhanced")
  .check_optional_pkg("ggnewscale", "ggraph_enhanced")

  dgi_edges$edge_type <- .normalise_interaction(dgi_edges$edge_type)

  gene_ids <- unique(dgi_edges$from)
  drug_ids <- unique(dgi_edges$to)

  drug_degree <- dplyr::count(dgi_edges, to, name = "targets")
  gene_degree <- dplyr::count(dgi_edges, from, name = "n_drugs")

  # ---- Build tidygraph ----
  nodes <- dplyr::bind_rows(
    tibble::tibble(name = gene_ids, type = "Gene") |>
      dplyr::left_join(gene_degree, by = c("name" = "from")) |>
      dplyr::mutate(degree = ifelse(is.na(n_drugs), 0L, n_drugs)) |>
      dplyr::select(name, type, degree),
    tibble::tibble(name = drug_ids, type = "Drug") |>
      dplyr::left_join(drug_degree, by = c("name" = "to")) |>
      dplyr::mutate(degree = ifelse(is.na(targets), 0L, targets)) |>
      dplyr::select(name, type, degree)
  ) |>
    dplyr::mutate(id = dplyr::row_number())

  edge_df <- dgi_edges |>
    dplyr::left_join(nodes, by = c("from" = "name")) |>
    dplyr::rename(from_id = id) |>
    dplyr::left_join(nodes, by = c("to" = "name")) |>
    dplyr::rename(to_id = id) |>
    dplyr::select(from = from_id, to = to_id, edge_type,
                  dplyr::any_of("score"))

  graph <- tidygraph::tbl_graph(nodes = nodes, edges = edge_df, directed = FALSE)

  # ---- Core mode: filter genes ----
  if (mode == "core") {
    core_gene_vec <- gene_degree |>
      dplyr::filter(n_drugs >= min_gene_edges) |>
      dplyr::pull(from)

    graph <- graph |>
      tidygraph::activate("nodes") |>
      dplyr::filter((type == "Gene" & name %in% core_gene_vec) | type == "Drug") |>
      tidygraph::activate("edges") |>
      dplyr::filter(!tidygraph::edge_is_multiple()) |>
      tidygraph::activate("nodes") |>
      dplyr::mutate(deg = tidygraph::centrality_degree()) |>
      dplyr::filter(deg > 0)
  }

  # ---- Compute qgraph FR layout and embed ----
  set.seed(seed)
  nd_tmp <- graph |> tidygraph::activate("nodes") |> tibble::as_tibble()
  ed_tmp <- graph |> tidygraph::activate("edges") |> tibble::as_tibble()
  # Map integer edge IDs back to node names
  ed_mapped <- ed_tmp |>
    dplyr::left_join(nd_tmp |> dplyr::select(id = .data$id, from_name = name),
                     by = c("from" = "id")) |>
    dplyr::left_join(nd_tmp |> dplyr::select(id = .data$id, to_name = name),
                     by = c("to" = "id")) |>
    dplyr::filter(!is.na(from_name) & !is.na(to_name)) |>
    dplyr::select(from = from_name, to = to_name)
  ig <- igraph::graph_from_data_frame(
    ed_mapped,
    vertices = nd_tmp |> dplyr::select(name),
    directed = FALSE
  )
  el <- igraph::as_edgelist(ig, names = FALSE)
  vc <- igraph::vcount(ig)
  fr <- qgraph::qgraph.layout.fruchtermanreingold(
    el, vcount = vc,
    area       = 2 * (vc^2.05),
    repulse.rad = (vc^3)
  )
  colnames(fr) <- c("x", "y")

  graph <- graph |>
    tidygraph::activate("nodes") |>
    dplyr::mutate(x = fr[, 1], y = fr[, 2])

  # ---- Node aesthetics ----
  graph <- graph |>
    tidygraph::activate("nodes") |>
    dplyr::mutate(
      is_hub     = (type == "Drug" & degree >= min_drug_targets) | (type == "Gene"),
      label_size = dplyr::if_else(type == "Gene", 3.2,
                           dplyr::if_else(degree >= 3, 2.5, 1.8)),
      label_face = dplyr::if_else(type == "Gene", "bold", "plain"),
      point_size = dplyr::case_when(
        type == "Gene" ~ scales_rescale(degree, to = c(3, 8)),
        type == "Drug" ~ scales_rescale(degree, to = c(1.5, 5))
      )
    )

  # ---- Edge palette ----
  edge_palette <- c(
    INHIBITOR  = "#E41A1C", AGONIST    = "#377EB8",
    ANTAGONIST = "#FF7F00", ACTIVATOR  = "#4DAF4A",
    BINDER     = "#984EA3", BLOCKER    = "#A65628",
    MODULATOR  = "#F781BF", SUBSTRATE  = "#66C2A5",
    None       = "#BBBBBB"
  )
  edge_types_present <- unique(dgi_edges$edge_type)
  edge_palette <- edge_palette[names(edge_palette) %in% edge_types_present]

  # ---- Plot ----
  nd <- graph |> tidygraph::activate("nodes") |> tibble::as_tibble()

  if (mode == "full") {
    label_df <- nd |> dplyr::filter(is_hub)
    subtitle_str <- paste0(
      sum(nd$type == "Gene"), " genes | ", sum(nd$type == "Drug"),
      " drugs | ", sum(nd$is_hub & nd$type == "Drug"),
      " hub drugs labelled | ", nrow(edge_df), " interactions"
    )
  } else {
    label_df <- nd
    subtitle_str <- paste0(
      "Genes with \u2265 ", min_gene_edges, " drug interactions | ",
      sum(nd$type == "Gene"), " genes, ", sum(nd$type == "Drug"),
      " drugs (all labelled)"
    )
  }

  p <- ggraph::ggraph(graph, layout = "manual", x = nd$x, y = nd$y) +
    ggraph::geom_edge_link(ggplot2::aes(color = edge_type),
                            alpha = 0.25, width = 0.3) +
    ggraph::scale_edge_colour_manual(
      values = edge_palette, name = "Interaction Type",
      guide  = ggplot2::guide_legend(override.aes = list(alpha = 0.8, width = 1.5))
    ) +
    ggnewscale::new_scale("colour") +
    ggraph::geom_node_point(ggplot2::aes(color = type, size = point_size),
                             alpha = 0.9) +
    ggplot2::scale_color_manual(
      values = c("Gene" = "#009E73", "Drug" = "#F0E442"),
      name = "Node Type"
    ) +
    ggplot2::scale_size_identity() +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(x = x, y = y, label = name,
                    size = I(label_size), fontface = label_face),
      max.overlaps = if (mode == "core") 200 else 100,
      family = "sans",
      segment.size = 0.2, segment.alpha = 0.5,
      box.padding = 0.3, force = if (mode == "core") 1 else 0.5,
      min.segment.length = if (mode == "core") 0.05 else 0.1
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle_str,
      caption  = paste0("DGIdb v5.0 | Fruchterman-Reingold layout | ", Sys.Date())
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = if (mode == "core") 20 else 18,
                                                face = "bold", hjust = 0.5),
      plot.subtitle    = ggplot2::element_text(size = if (mode == "core") 11 else 10,
                                                hjust = 0.5, color = "grey40"),
      plot.caption     = ggplot2::element_text(size = 8, color = "grey60"),
      legend.position  = "bottom",
      legend.box       = "vertical",
      legend.text      = ggplot2::element_text(size = 9),
      legend.title     = ggplot2::element_text(size = 10, face = "bold"),
      legend.key.size  = ggplot2::unit(0.4, "cm"),
      plot.margin      = ggplot2::margin(15, 15, 15, 15),
      plot.background  = ggplot2::element_rect(fill = "white", color = NA)
    )

  .save_ggplot(p, out_dir, base_filename, width, height, dpi)
  invisible(p)
}
