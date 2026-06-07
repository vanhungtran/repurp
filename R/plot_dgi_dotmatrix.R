#' Plot a Gene × Drug-Class dot matrix
#'
#' Renders a publication-quality dot matrix (also called a balloon plot or
#' bubble heatmap) where rows are genes, columns are drug classes, dot size
#' encodes the number of drug–gene connections within that class, and dot colour
#' encodes the dominant interaction type. This plot gives an at-a-glance
#' pharmacological summary useful for drug repurposing.
#'
#' @param dgi_edges A data frame with columns `from` (gene), `to` (drug),
#'   `edge_type` (interaction type string).
#' @param drug_info A data frame with columns `drug_name` and `category`.
#'   Defaults to [repurp_ad_drugs()].
#' @param title Plot title.
#' @param max_genes Maximum number of genes shown (top by total degree).
#'   Default `40L`.
#' @param min_connections Minimum total drug–gene connections for a cell to be
#'   drawn. Default `1L`.
#' @param edge_palette Named hex vector for interaction types. Defaults to
#'   [dgi_edge_colors].
#' @param out_dir Directory for saved files. `NULL` = plot to current device only.
#' @param base_filename Base filename without extension.
#' @param width,height Plot dimensions in inches. Defaults `14` / `12`.
#' @param dpi Resolution for PNG. Default `300`.
#'
#' @return Invisibly returns the `ggplot` object.
#' @export
#' @examples
#' \dontrun{
#' dgi <- dplyr::transmute(net_dat, from = Gene, to = Drug, edge_type = InteractionType)
#' plot_dgi_dotmatrix(dgi, drug_info = repurp_ad_drugs(), out_dir = "results/")
#' }
plot_dgi_dotmatrix <- function(dgi_edges,
                                drug_info     = repurp_ad_drugs(),
                                title         = "Gene × Drug-Class Connectivity Matrix",
                                max_genes     = 40L,
                                min_connections = 1L,
                                edge_palette  = dgi_edge_colors,
                                out_dir       = NULL,
                                base_filename = paste0("DGI_dotmatrix_", Sys.Date()),
                                width  = 14,
                                height = 12,
                                dpi    = 300) {

  .check_optional_pkg("ggplot2", "plot_dgi_dotmatrix")

  dgi_edges$edge_type <- .normalise_interaction(dgi_edges$edge_type)

  # Attach drug class
  cat_lookup <- stats::setNames(drug_info$category, drug_info$drug_name)
  dgi_edges$drug_class <- cat_lookup[dgi_edges$to]
  dgi_edges$drug_class[is.na(dgi_edges$drug_class)] <- "Other"

  # Top genes by total degree
  top_genes <- dgi_edges |>
    dplyr::count(from, sort = TRUE) |>
    dplyr::slice_head(n = max_genes) |>
    dplyr::pull(from)

  df <- dgi_edges |>
    dplyr::filter(from %in% top_genes)

  # Aggregate: count connections and dominant interaction type per gene × drug_class
  cell_df <- df |>
    dplyr::count(gene = from, drug_class, edge_type, name = "n") |>
    dplyr::group_by(gene, drug_class) |>
    dplyr::mutate(total = sum(n)) |>
    dplyr::slice_max(n, n = 1L, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::filter(total >= min_connections)

  # Ordered axes
  gene_order  <- cell_df |>
    dplyr::count(gene, wt = total, sort = TRUE) |>
    dplyr::pull(gene)
  class_order <- cell_df |>
    dplyr::count(drug_class, wt = total, sort = TRUE) |>
    dplyr::pull(drug_class)

  cell_df$gene       <- factor(cell_df$gene,       levels = rev(gene_order))
  cell_df$drug_class <- factor(cell_df$drug_class, levels = class_order)

  # Colour palette
  ep <- edge_palette[names(edge_palette) %in% cell_df$edge_type]
  missing_types <- setdiff(unique(cell_df$edge_type), names(ep))
  if (length(missing_types) > 0) {
    ep <- c(ep, stats::setNames(rep("#BBBBBB", length(missing_types)), missing_types))
  }

  max_total <- max(cell_df$total, na.rm = TRUE)

  p <- ggplot2::ggplot(cell_df,
    ggplot2::aes(x = drug_class, y = gene,
                 size = total, color = edge_type)) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_size_area(
      max_size = 12,
      name     = "# connections",
      breaks   = scales_rescale(
        c(1, 0.33, 0.66, 1) * max_total,
        to = c(1, max_total)
      ) |> round() |> unique()
    ) +
    ggplot2::scale_color_manual(
      values = ep,
      name   = "Dominant interaction",
      guide  = ggplot2::guide_legend(
        override.aes = list(size = 4),
        ncol = 1
      )
    ) +
    # Count label inside dot
    ggplot2::geom_text(
      ggplot2::aes(label = total),
      size = 2.5, color = "white", fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = paste0(
        length(unique(cell_df$gene)), " genes × ",
        length(unique(cell_df$drug_class)), " drug classes | ",
        "size = # drug–gene connections | colour = dominant interaction type"
      ),
      caption = paste0("DGIdb v5.0 | Repurp package | ", Sys.Date()),
      x = "Drug class", y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle   = ggplot2::element_text(size = 8,  hjust = 0.5, color = "grey40"),
      plot.caption    = ggplot2::element_text(size = 7,  color = "grey60"),
      axis.text.x     = ggplot2::element_text(size = 9,  angle = 35, hjust = 1,
                                               face = "bold"),
      axis.text.y     = ggplot2::element_text(size = 8,  color = "#009E73",
                                               face = "bold"),
      panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "right",
      legend.text      = ggplot2::element_text(size = 8),
      legend.title     = ggplot2::element_text(size = 9, face = "bold"),
      plot.background  = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin      = ggplot2::margin(10, 10, 10, 10)
    )

  .save_ggplot(p, out_dir, base_filename, width, height, dpi)
  invisible(p)
}
