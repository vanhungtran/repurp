#' Plot a circular arc / edge-bundle drug–gene network
#'
#' Draws a circular network where drugs and genes alternate around the
#' circumference, connected by curved arc edges. Drug edges are coloured by
#' drug class. This layout is compact and well-suited for posters or
#' presentations where space is tight.
#'
#' @param edges A data frame with columns `from` (drug), `to` (gene), and
#'   `class` (drug class string, used for edge colouring). Typically built by
#'   joining [repurp_ad_edges()] with [repurp_ad_drugs()].
#' @param title Plot title.
#' @param drug_class_colors Named character vector mapping drug class strings
#'   to hex colours. Defaults to [repurp_class_colors()].
#' @param gene_color Hex colour for gene nodes. Default `"#009E73"`.
#' @param drug_color Hex colour for drug nodes. Default `"#E41A1C"`.
#' @param seed Random seed for reproducible layout. Default `42`.
#' @param out_dir Directory for saved files. `NULL` = draw to current device.
#' @param base_filename Base filename without extension.
#' @param width,height Plot dimensions in inches. Defaults `16` / `16`.
#' @param dpi Resolution for PNG. Default `300`.
#'
#' @return Invisibly returns the `ggplot` object.
#' @export
#' @examples
#' \dontrun{
#' edges <- repurp_ad_edges() |>
#'   dplyr::filter(to %in% my_genes) |>
#'   dplyr::left_join(repurp_ad_drugs(), by = c("from" = "drug_name")) |>
#'   dplyr::rename(class = category)
#'
#' plot_dgi_circular_arc(edges, out_dir = "results/")
#' }
plot_dgi_circular_arc <- function(edges,
                                   title             = "Drug\u2013Biomarker Interaction Network",
                                   drug_class_colors = repurp_class_colors(),
                                   gene_color        = "#009E73",
                                   drug_color        = "#E41A1C",
                                   seed              = 42L,
                                   out_dir           = NULL,
                                   base_filename     = paste0("DGI_circular_arc_", Sys.Date()),
                                   width  = 16,
                                   height = 16,
                                   dpi    = 300) {

  .check_ggraph()

  # ---- Build node table ----
  drug_ids <- unique(edges$from)
  gene_ids <- unique(edges$to)

  nodes <- dplyr::bind_rows(
    tibble::tibble(name = drug_ids, type = "Drug"),
    tibble::tibble(name = gene_ids, type = "Biomarker")
  ) |>
    dplyr::mutate(id = dplyr::row_number())

  # ---- Build igraph from edges ----
  ig <- igraph::graph_from_data_frame(
    edges |> dplyr::select(from, to, class),
    vertices = nodes |> dplyr::select(name),
    directed = FALSE
  )

  graph_tbl <- tidygraph::as_tbl_graph(ig)

  # ---- Edge colour palette ----
  edge_pal <- drug_class_colors
  edge_pal <- edge_pal[names(edge_pal) %in% unique(edges$class)]

  # ---- Plot ----
  set.seed(seed)
  p <- ggraph::ggraph(graph_tbl, layout = "linear", circular = TRUE) +
    ggraph::geom_edge_arc(
      ggplot2::aes(color = class),
      alpha = 0.25, strength = 0.2, width = 0.5
    ) +
    ggraph::scale_edge_colour_manual(
      values = edge_pal,
      guide = "none"
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(color = type, size = type),
      alpha = 0.9
    ) +
    ggplot2::scale_color_manual(
      values = c("Drug" = drug_color, "Biomarker" = gene_color),
      name = "Node Type"
    ) +
    ggplot2::scale_size_manual(
      values = c("Drug" = 4, "Biomarker" = 2.5),
      guide = "none"
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(
        x = x * 1.05,
        y = y * 1.05,
        label = name,
        color = type,
        angle = dplyr::if_else(
          atan2(y, x) * 180 / pi < -90 | atan2(y, x) * 180 / pi > 90,
          atan2(y, x) * 180 / pi + 180,
          atan2(y, x) * 180 / pi
        ),
        hjust = dplyr::if_else(
          atan2(y, x) * 180 / pi < -90 | atan2(y, x) * 180 / pi > 90,
          1, 0
        )
      ),
      size = 3.5, fontface = "bold"
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = paste0(
        length(drug_ids), " drugs | ",
        length(gene_ids), " biomarkers | ",
        nrow(edges), " curated edges"
      ),
      caption  = paste0("Curated mechanism-of-action | Circular arc layout | ",
                        Sys.Date())
    ) +
    ggplot2::scale_x_continuous(limits = c(-2.5, 2.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(-2.5, 2.5), expand = c(0, 0)) +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 18, face = "bold",
                                                hjust = 0.5),
      plot.subtitle    = ggplot2::element_text(size = 10, hjust = 0.5,
                                                color = "grey40"),
      plot.caption     = ggplot2::element_text(size = 7, color = "grey60"),
      legend.position  = "bottom",
      plot.margin      = ggplot2::margin(1.5, 1.5, 1.5, 1.5, "in"),
      plot.background  = ggplot2::element_rect(fill = "white", color = NA)
    )

  .save_ggplot(p, out_dir, base_filename, width, height, dpi)
  invisible(p)
}
