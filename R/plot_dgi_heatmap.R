#' Plot a Drug × Gene interaction tile heatmap
#'
#' Renders a publication-quality ggplot2 tile heatmap where rows are drugs,
#' columns are genes, and fill encodes the dominant interaction type.
#' An optional left-side strip encodes drug class when `drug_info` is supplied.
#' Rows are sorted by drug class then degree; columns by gene degree.
#'
#' @param dgi_edges A data frame with columns `from` (gene), `to` (drug),
#'   `edge_type` (interaction type string).
#' @param drug_info Optional data frame with columns `drug_name` and `category`.
#'   When supplied a colour strip on the left side shows drug class.
#'   Defaults to [repurp_ad_drugs()] if all drugs in `dgi_edges` match AD drugs;
#'   set to `NULL` to suppress the strip entirely.
#' @param title Plot title.
#' @param max_drugs Maximum number of drugs shown (top by degree). Default `60L`.
#' @param max_genes Maximum number of genes shown (top by degree). Default `50L`.
#' @param edge_palette Named character vector of interaction type → hex colour.
#'   Defaults to [dgi_edge_colors].
#' @param out_dir Directory for saved files. `NULL` = plot to current device only.
#' @param base_filename Base filename without extension.
#' @param width,height Plot dimensions in inches. Defaults `16` / `14`.
#' @param dpi Resolution for PNG. Default `300`.
#'
#' @return Invisibly returns the `ggplot` object.
#' @export
#' @examples
#' \dontrun{
#' dgi <- dplyr::transmute(net_dat, from = Gene, to = Drug, edge_type = InteractionType)
#' plot_dgi_heatmap(dgi, drug_info = repurp_ad_drugs(), out_dir = "results/")
#' }
plot_dgi_heatmap <- function(dgi_edges,
                              drug_info      = NULL,
                              title          = "Drug–Gene Interaction Heatmap",
                              max_drugs      = 60L,
                              max_genes      = 50L,
                              edge_palette   = dgi_edge_colors,
                              out_dir        = NULL,
                              base_filename  = paste0("DGI_heatmap_", Sys.Date()),
                              width  = 16,
                              height = 14,
                              dpi    = 300) {

  .check_optional_pkg("ggplot2", "plot_dgi_heatmap")

  dgi_edges$edge_type <- .normalise_interaction(dgi_edges$edge_type)

  # Select top drugs and genes by degree
  top_drugs <- dgi_edges |>
    dplyr::count(to, sort = TRUE) |>
    dplyr::slice_head(n = max_drugs) |>
    dplyr::pull(to)

  top_genes <- dgi_edges |>
    dplyr::count(from, sort = TRUE) |>
    dplyr::slice_head(n = max_genes) |>
    dplyr::pull(from)

  df <- dgi_edges |>
    dplyr::filter(to %in% top_drugs, from %in% top_genes)

  # Dominant interaction type per drug-gene pair
  tile_df <- df |>
    dplyr::count(drug = to, gene = from, edge_type, sort = TRUE) |>
    dplyr::group_by(drug, gene) |>
    dplyr::slice_max(n, n = 1L, with_ties = FALSE) |>
    dplyr::ungroup()

  # Attach drug category for sorting / strip
  if (!is.null(drug_info)) {
    cat_lookup <- stats::setNames(drug_info$category, drug_info$drug_name)
    tile_df$category <- cat_lookup[tile_df$drug]
    tile_df$category[is.na(tile_df$category)] <- "Other"
  } else {
    tile_df$category <- "Drug"
  }

  # Drug order: by category then degree
  drug_order <- tile_df |>
    dplyr::count(drug, category, sort = FALSE) |>
    dplyr::arrange(category, dplyr::desc(n)) |>
    dplyr::pull(drug)

  # Gene order: by degree desc
  gene_order <- tile_df |>
    dplyr::count(gene, sort = TRUE) |>
    dplyr::pull(gene)

  tile_df$drug <- factor(tile_df$drug, levels = rev(drug_order))
  tile_df$gene <- factor(tile_df$gene, levels = gene_order)

  # Colour palette — keep only types present
  ep <- edge_palette[names(edge_palette) %in% tile_df$edge_type]
  missing_types <- setdiff(unique(tile_df$edge_type), names(ep))
  if (length(missing_types) > 0) {
    extras <- stats::setNames(rep("#BBBBBB", length(missing_types)), missing_types)
    ep <- c(ep, extras)
  }

  # Optional category strip colours (left side via geom_tile in a secondary layer)
  cat_levels <- if (!is.null(drug_info)) {
    unique(tile_df$category[order(tile_df$drug)])
  } else "Drug"

  # --- Build plot ---
  p <- ggplot2::ggplot(tile_df,
    ggplot2::aes(x = gene, y = drug, fill = edge_type)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(
      values = ep,
      name   = "Interaction type",
      guide  = ggplot2::guide_legend(
        ncol = 2,
        override.aes = list(color = NA)
      )
    )

  # Category strip: thin coloured rect on the left
  if (!is.null(drug_info)) {
    strip_df <- tile_df |>
      dplyr::distinct(drug, category) |>
      dplyr::mutate(x_strip = -0.9)

    # Use ad_drug_class_colors if available, else generate palette
    n_cats <- length(unique(strip_df$category))
    cat_pal <- grDevices::hcl.colors(n_cats, palette = "Set2")
    cat_pal <- stats::setNames(cat_pal, unique(strip_df$category))

    p <- p +
      ggplot2::geom_tile(
        data = strip_df,
        ggplot2::aes(x = x_strip, y = drug, fill = NULL, color = NULL),
        fill  = cat_pal[strip_df$category],
        width = 0.7, color = NA
      ) +
      ggplot2::annotate("text",
        x = -0.9,
        y = length(levels(tile_df$drug)) / 2 + 0.5,
        label = "", size = 0
      )
  }

  n_drugs <- length(levels(tile_df$drug))
  n_genes <- length(levels(tile_df$gene))

  p <- p +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = c(1.5, 0.5))) +
    ggplot2::labs(
      title    = title,
      subtitle = paste0(n_drugs, " drugs × ", n_genes, " genes | ",
                        "fill = dominant interaction type"),
      caption  = paste0("DGIdb v5.0 | Repurp package | ", Sys.Date()),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle    = ggplot2::element_text(size = 9, hjust = 0.5, color = "grey40"),
      plot.caption     = ggplot2::element_text(size = 7, color = "grey60"),
      axis.text.x      = ggplot2::element_text(size = 7, angle = 45, hjust = 1,
                                                color = "#009E73", face = "bold"),
      axis.text.y      = ggplot2::element_text(size = 7),
      panel.grid       = ggplot2::element_blank(),
      legend.position  = "right",
      legend.text      = ggplot2::element_text(size = 8),
      legend.title     = ggplot2::element_text(size = 9, face = "bold"),
      plot.background  = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin      = ggplot2::margin(10, 10, 10, 10)
    )

  .save_ggplot(p, out_dir, base_filename, width, height, dpi)
  invisible(p)
}
