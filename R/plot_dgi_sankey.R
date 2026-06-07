#' Plot a Drug Class → Drug → Gene alluvial (Sankey) flow diagram
#'
#' Renders a three-axis alluvial plot showing how drug classes funnel through
#' individual drugs into biomarker genes. The width of each flow band is
#' proportional to the number of drug–gene connections. Requires the
#' `ggalluvial` and `ggplot2` packages.
#'
#' @param dgi_edges A data frame with columns `from` (gene), `to` (drug),
#'   `edge_type` (interaction type string).
#' @param drug_info A data frame with columns `drug_name` and `category`.
#'   Defaults to [repurp_ad_drugs()].
#' @param title Plot title.
#' @param min_connections Minimum number of edges for a drug–gene pair to be
#'   included. Increase to reduce clutter in large networks. Default `1L`.
#' @param max_genes Maximum number of genes shown (top by degree). Default `30L`.
#' @param fill_by One of `"drug_class"` (colour alluvium by drug category) or
#'   `"edge_type"` (colour by interaction type). Default `"drug_class"`.
#' @param edge_palette Named hex vector for interaction types (used when
#'   `fill_by = "edge_type"`). Defaults to [dgi_edge_colors].
#' @param out_dir Directory for saved files. `NULL` = plot to current device only.
#' @param base_filename Base filename without extension.
#' @param width,height Plot dimensions in inches. Defaults `18` / `14`.
#' @param dpi Resolution for PNG. Default `300`.
#'
#' @return Invisibly returns the `ggplot` object.
#' @export
#' @examples
#' \dontrun{
#' dgi <- dplyr::transmute(net_dat, from = Gene, to = Drug, edge_type = InteractionType)
#' plot_dgi_sankey(dgi, drug_info = repurp_ad_drugs(), out_dir = "results/")
#' }
plot_dgi_sankey <- function(dgi_edges,
                             drug_info       = repurp_ad_drugs(),
                             title           = "Drug Class → Drug → Gene Flow",
                             min_connections = 1L,
                             max_genes       = 30L,
                             fill_by         = c("drug_class", "edge_type"),
                             edge_palette    = dgi_edge_colors,
                             out_dir         = NULL,
                             base_filename   = paste0("DGI_sankey_", Sys.Date()),
                             width  = 18,
                             height = 14,
                             dpi    = 300) {

  fill_by <- match.arg(fill_by)
  .check_optional_pkg("ggplot2",   "plot_dgi_sankey")
  .check_optional_pkg("ggalluvial","plot_dgi_sankey")

  dgi_edges$edge_type <- .normalise_interaction(dgi_edges$edge_type)

  # Attach drug class
  cat_lookup <- stats::setNames(drug_info$category, drug_info$drug_name)
  dgi_edges$drug_class <- cat_lookup[dgi_edges$to]
  dgi_edges$drug_class[is.na(dgi_edges$drug_class)] <- "Other"

  # Limit to top genes by degree
  top_genes <- dgi_edges |>
    dplyr::count(from, sort = TRUE) |>
    dplyr::slice_head(n = max_genes) |>
    dplyr::pull(from)

  df <- dgi_edges |>
    dplyr::filter(from %in% top_genes) |>
    dplyr::count(drug_class, drug = to, gene = from, edge_type,
                 name = "freq") |>
    dplyr::filter(freq >= min_connections)

  # Ordered axes so that many-edge items sit centrally
  class_order <- df |>
    dplyr::count(drug_class, wt = freq, sort = TRUE) |>
    dplyr::pull(drug_class)
  drug_order  <- df |>
    dplyr::count(drug, wt = freq, sort = TRUE) |>
    dplyr::pull(drug)
  gene_order  <- df |>
    dplyr::count(gene, wt = freq, sort = TRUE) |>
    dplyr::pull(gene)

  df$drug_class <- factor(df$drug_class, levels = class_order)
  df$drug       <- factor(df$drug,       levels = drug_order)
  df$gene       <- factor(df$gene,       levels = gene_order)

  # Colour scale
  if (fill_by == "drug_class") {
    n_cls   <- length(class_order)
    pal     <- grDevices::hcl.colors(n_cls, palette = "Dark 3")
    fill_col <- stats::setNames(pal, class_order)
    fill_var <- "drug_class"
    legend_title <- "Drug class"
  } else {
    ep <- edge_palette[names(edge_palette) %in% df$edge_type]
    fill_col <- ep
    fill_var <- "edge_type"
    legend_title <- "Interaction type"
  }

  p <- ggplot2::ggplot(df,
    ggplot2::aes(
      axis1 = drug_class,
      axis2 = drug,
      axis3 = gene,
      y     = freq
    )
  ) +
    ggalluvial::geom_alluvium(
      ggplot2::aes(fill = .data[[fill_var]]),
      alpha = 0.65, width = 1/3, knot.pos = 0.4
    ) +
    ggalluvial::geom_stratum(
      fill  = "grey92", color = "grey60",
      width = 1/3, linewidth = 0.3
    ) +
    ggalluvial::geom_flow(
      ggplot2::aes(fill = .data[[fill_var]]),
      alpha = 0.4, width = 1/3
    ) +
    ggplot2::geom_text(
      stat = ggalluvial::StatStratum,
      ggplot2::aes(label = ggplot2::after_stat(stratum)),
      size = 2.8, fontface = "bold", color = "grey20"
    ) +
    ggplot2::scale_fill_manual(
      values = fill_col,
      name   = legend_title,
      guide  = ggplot2::guide_legend(ncol = 2)
    ) +
    ggplot2::scale_x_discrete(
      limits = c("Drug Class", "Drug", "Gene"),
      expand = ggplot2::expansion(mult = c(0.05, 0.05))
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = paste0(
        length(unique(df$drug_class)), " drug classes | ",
        length(unique(df$drug)),       " drugs | ",
        length(unique(df$gene)),        " genes"
      ),
      caption = paste0("DGIdb v5.0 | Repurp package | ", Sys.Date()),
      x = NULL, y = "Connections"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle   = ggplot2::element_text(size = 9,  hjust = 0.5, color = "grey40"),
      plot.caption    = ggplot2::element_text(size = 7,  color = "grey60"),
      axis.text.x     = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y     = ggplot2::element_blank(),
      axis.ticks.y    = ggplot2::element_blank(),
      panel.grid      = ggplot2::element_blank(),
      legend.position = "right",
      legend.text     = ggplot2::element_text(size = 8),
      legend.title    = ggplot2::element_text(size = 9, face = "bold"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin     = ggplot2::margin(15, 15, 15, 15)
    )

  .save_ggplot(p, out_dir, base_filename, width, height, dpi)
  invisible(p)
}
