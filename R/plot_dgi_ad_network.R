#' Plot an AD treatment drug–biomarker network with drug categories
#'
#' Builds a publication-quality network showing established atopic-dermatitis (AD)
#' treatments (source: [repurp_ad_drugs()]) connected to biomarker genes via
#' curated mechanism-of-action edges (source: [repurp_ad_edges()]).
#'
#' The plot uses `igraph` + `qgraph` Fruchterman-Reingold layout, rendered with
#' `ggplot2`, `ggrepel` labels, and `ggnewscale` for separate legends.
#' Drug nodes are coloured by therapeutic category; biomarker nodes are green.
#' Solid edges = direct molecular target / DGIdb confirmed;
#' dashed edges = pathway-level (indirect) regulation.
#'
#' @param edges A data frame with columns `from` (drug), `to` (gene),
#'   `edge_type` (`"direct"`, `"indirect"`, or a DGIdb interaction type string),
#'   and `mechanism` (description string). Typically built by combining
#'   [repurp_ad_edges()] with optional DGIdb matches.
#' @param drug_info A data frame with columns `drug_name`, `brand`, `category`.
#'   Defaults to [repurp_ad_drugs()].
#' @param title Plot title.
#' @param gene_color Hex colour for biomarker gene nodes. Default `"#009E73"`.
#' @param out_dir Directory for saved files. `NULL` = draw to current device only.
#' @param base_filename Base filename without extension.
#' @param width,height Plot dimensions in inches. Defaults `18` / `16`.
#' @param dpi Resolution for PNG. Default `300`.
#'
#' @return Invisibly returns the `ggplot` object.
#' @export
#' @examples
#' \dontrun{
#' edges <- repurp_ad_edges() |>
#'   dplyr::filter(to %in% my_genes)
#'
#' plot_dgi_ad_network(
#'   edges    = edges,
#'   drug_info = repurp_ad_drugs(),
#'   out_dir  = "results/"
#' )
#' }
plot_dgi_ad_network <- function(edges,
                                 drug_info      = repurp_ad_drugs(),
                                 title          = "AD Treatment Drug\u2013Biomarker Network",
                                 gene_color     = "#009E73",
                                 out_dir        = NULL,
                                 base_filename  = paste0("AD_drug_network_", Sys.Date()),
                                 width  = 18,
                                 height = 16,
                                 dpi    = 300) {

  .check_optional_pkg("ggraph", "ad_network")
  .check_optional_pkg("tidygraph", "ad_network")
  .check_optional_pkg("ggrepel", "ad_network")
  .check_optional_pkg("ggnewscale", "ad_network")

  # ---- Build node table ----
  drug_names <- unique(edges$from)
  gene_names <- unique(edges$to)

  drug_cat <- drug_info |>
    dplyr::select(name = drug_name, category, brand) |>
    dplyr::filter(name %in% drug_names)

  nodes <- dplyr::bind_rows(
    tibble::tibble(
      name   = drug_names,
      type   = "AD Drug",
      group  = drug_cat$category[match(drug_names, drug_cat$name)],
      label2 = drug_cat$brand[match(drug_names, drug_cat$name)]
    ),
    tibble::tibble(
      name   = gene_names,
      type   = "Biomarker",
      group  = "Biomarker",
      label2 = ""
    )
  ) |>
    dplyr::mutate(id = dplyr::row_number())

  # ---- Edge table with node IDs ----
  edge_df <- edges |>
    dplyr::left_join(nodes, by = c("from" = "name")) |>
    dplyr::rename(from_id = id) |>
    dplyr::left_join(nodes, by = c("to" = "name")) |>
    dplyr::rename(to_id = id) |>
    dplyr::mutate(
      edge_color = dplyr::case_when(
        grepl("DGIdb", edge_type)         ~ "#E41A1C",
        edge_type == "direct"             ~ "#377EB8",
        TRUE                              ~ "#999999"
      ),
      edge_linetype = dplyr::if_else(
        grepl("DGIdb|direct", edge_type), "solid", "dashed"
      ),
      edge_alpha = dplyr::if_else(
        grepl("DGIdb|direct", edge_type), 0.8, 0.4
      ),
      edge_width = dplyr::if_else(
        grepl("DGIdb|direct", edge_type), 0.8, 0.4
      )
    )

  # ---- Build graph for layout computation ----
  graph <- tidygraph::tbl_graph(
    nodes = nodes,
    edges = edge_df |> dplyr::select(from = from_id, to = to_id,
                                      dplyr::any_of(c("edge_type", "mechanism",
                                                       "edge_color", "edge_linetype",
                                                       "edge_alpha", "edge_width"))),
    directed = FALSE
  )

  # Compute qgraph FR layout
  nd_tmp <- graph |> tidygraph::activate("nodes") |> tibble::as_tibble()
  ed_tmp <- graph |> tidygraph::activate("edges") |> tibble::as_tibble()
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
    area       = 12 * (vc^2),
    repulse.rad = (vc^2.8)
  )

  nodes$x <- fr[, 1]
  nodes$y <- fr[, 2]

  # ---- Node aesthetics ----
  drug_degree <- edges |>
    dplyr::count(from, name = "n_connections")

  cat_colors <- c(
    "Corticosteroid"           = "#FDE0DD",
    "Mild topical steroid"     = "#FDE0DD",
    "Medium topical steroid"   = "#F4B6C2",
    "Strong topical steroid"   = "#E55970",
    "Topical steroid"          = "#E41A1C",
    "Calcineurin inhibitor"    = "#FDD49E",
    "PDE4 inhibitor"           = "#AED6F1",
    "JAK inhibitor"            = "#A9DFBF",
    "JAK inhibitor (topical)"  = "#A9DFBF",
    "JAK inhibitor (oral)"     = "#7DCEA0",
    "Biologic"                 = "#D7BDE2",
    "Biologic (anti-IL4R)"     = "#D7BDE2",
    "Biologic (anti-IL13)"     = "#BB8FCE",
    "Biologic (anti-IL31RA)"   = "#A569BD",
    "AhR agonist"              = "#F9E79F",
    "Immunosuppressant"        = "#F5B7B1"
  )

  nodes <- nodes |>
    dplyr::left_join(drug_degree, by = c("name" = "from")) |>
    dplyr::mutate(
      n_connections = ifelse(is.na(n_connections), 1, n_connections),
      point_size = dplyr::case_when(
        type == "AD Drug"   ~ scales_rescale(n_connections, to = c(6, 14)),
        type == "Biomarker" ~ scales_rescale(n_connections, to = c(3, 8))
      ),
      label = dplyr::if_else(type == "AD Drug" & !is.na(label2) & label2 != "",
                              paste0(name, "\n(", label2, ")"), name),
      label_size = dplyr::if_else(type == "AD Drug", 3.0, 2.5),
      label_face = dplyr::if_else(type == "AD Drug", "bold", "plain"),
      node_fill = dplyr::case_when(
        type == "AD Drug"   ~ cat_colors[group],
        type == "Biomarker" ~ gene_color
      )
    )

  # ---- Build edge coordinate data ----
  edge_coords <- edge_df |>
    dplyr::left_join(nodes |> dplyr::select(name, x_from = x, y_from = y),
                      by = c("from" = "name")) |>
    dplyr::left_join(nodes |> dplyr::select(name, x_to = x, y_to = y),
                      by = c("to" = "name"))

  # ---- Legend data ----
  legend_df <- nodes |>
    dplyr::filter(type == "AD Drug") |>
    dplyr::distinct(name, group, node_fill)

  # ---- Plot ----
  p <- ggplot2::ggplot() +
    # Edges
    ggplot2::geom_segment(
      data = edge_coords,
      ggplot2::aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
                    color = edge_type, linetype = I(edge_linetype),
                    linewidth = edge_width, alpha = I(edge_alpha))
    ) +
    # Biomarker nodes
    ggplot2::geom_point(
      data = nodes |> dplyr::filter(type == "Biomarker"),
      ggplot2::aes(x = x, y = y, size = point_size),
      color = gene_color, alpha = 0.9
    ) +
    # Drug nodes (filled squares with dark border)
    ggplot2::geom_point(
      data = nodes |> dplyr::filter(type == "AD Drug"),
      ggplot2::aes(x = x, y = y, size = point_size, fill = node_fill),
      color = "grey30", stroke = 1.2, shape = 22, alpha = 0.95
    ) +
    ggplot2::scale_fill_identity() +
    # Biomarker labels
    ggrepel::geom_text_repel(
      data = nodes |> dplyr::filter(type == "Biomarker"),
      ggplot2::aes(x = x, y = y, label = name, size = I(label_size),
                    fontface = label_face),
      max.overlaps = 50, family = "sans", color = gene_color,
      segment.size = 0.2, segment.alpha = 0.4,
      box.padding = 0.3, force = 1, min.segment.length = 0.1
    ) +
    # Drug labels
    ggrepel::geom_text_repel(
      data = nodes |> dplyr::filter(type == "AD Drug"),
      ggplot2::aes(x = x, y = y, label = label, size = I(label_size),
                    fontface = label_face),
      max.overlaps = 50, family = "sans", color = "grey20",
      segment.size = 0.2, segment.alpha = 0.6,
      box.padding = 0.4, force = 2, min.segment.length = 0.1
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_color_manual(
      values = c(
        direct   = "#377EB8",
        indirect = "#999999",
        "DGIdb: INHIBITOR" = "#E41A1C",
        "DGIdb: AGONIST"   = "#377EB8",
        "DGIdb: ACTIVATOR" = "#4DAF4A",
        "DGIdb: BINDER"    = "#984EA3",
        "DGIdb: None"      = "#999999"
      ),
      name = "Connection type"
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = paste0(
        length(drug_names), " AD drugs (colored by class) | ",
        length(gene_names), " connected biomarkers | ",
        nrow(edges), " edges"
      ),
      caption  = paste0(
        "Solid edges = DGIdb confirmed or direct molecular target | ",
        "Dashed = pathway-level regulation | ", Sys.Date()
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 18, face = "bold",
                                                hjust = 0.5),
      plot.subtitle    = ggplot2::element_text(size = 9, hjust = 0.5,
                                                color = "grey40"),
      plot.caption     = ggplot2::element_text(size = 7, color = "grey60"),
      legend.position  = "bottom",
      legend.text      = ggplot2::element_text(size = 8),
      legend.title     = ggplot2::element_text(size = 9, face = "bold"),
      legend.key.size  = ggplot2::unit(0.4, "cm"),
      plot.margin      = ggplot2::margin(15, 15, 15, 15),
      plot.background  = ggplot2::element_rect(fill = "white", color = NA)
    )

  .save_ggplot(p, out_dir, base_filename, width, height, dpi)
  invisible(p)
}
