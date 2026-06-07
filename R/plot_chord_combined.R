#' Draw a combined chord diagram of curated AD drugs + DGIdb hub drugs
#'
#' Builds a [circlize](https://cran.r-project.org/package=circlize) chord diagram
#' that merges (1) curated AD drug–biomarker edges from [repurp_ad_edges()] and
#' (2) DGIdb non-AD drugs filtered to those with at least `min_dgidb_targets`
#' gene connections. Drug sectors are coloured by drug class / interaction type;
#' gene sectors are green. Link colour encodes interaction type.
#'
#' @inheritParams plot_chord_dgi
#' @param dgi_edges A data frame with columns `drug`, `gene`, `int_type`
#'   (interaction type). Typically built from DGIdb non-AD drug data.
#' @param curated_edges A data frame with columns `drug`, `gene`, `class`
#'   (drug class), and `int_type` (set to `"Curated MOA"`). Typically from
#'   [repurp_ad_edges()] joined with [repurp_ad_drugs()].
#' @param min_dgidb_targets Minimum number of gene targets for a DGIdb drug to be
#'   included. Default `2`.
#'
#' @return Invisibly returns the adjacency matrix used for plotting.
#' @export
#' @examples
#' \dontrun{
#' curated <- repurp_ad_edges() |>
#'   dplyr::filter(to %in% my_genes) |>
#'   dplyr::left_join(repurp_ad_drugs(), by = c("from" = "drug_name")) |>
#'   dplyr::rename(drug = from, gene = to, class = category) |>
#'   dplyr::mutate(int_type = "Curated MOA")
#'
#' dgi <- dplyr::transmute(my_dgidb_data,
#'   drug = Drug, gene = Gene, int_type = InteractionType
#' )
#'
#' plot_chord_dgi_combined(curated_edges = curated, dgi_edges = dgi, out_dir = "results/")
#' }
plot_chord_dgi_combined <- function(dgi_edges,
                                     curated_edges,
                                     drug_class_colors = repurp_class_colors(),
                                     gene_color        = "#009E73",
                                     title             = "Drug–Gene Chord Diagram — AD + DGIdb",
                                     out_dir           = NULL,
                                     base_filename     = paste0("DGI_chord_combined_", Sys.Date()),
                                     svg_w             = 32,
                                     svg_h             = 32,
                                     transparency      = 0.15,
                                     min_dgidb_targets = 2L) {

  if (!requireNamespace("circlize", quietly = TRUE)) {
    rlang::abort("Package 'circlize' is required. Install with: install.packages('circlize')")
  }

  # ---- Filter DGIdb drugs by hub status ----
  dgi_edges$int_type <- .normalise_interaction(dgi_edges$int_type)
  dgi_edges <- dgi_edges |>
    dplyr::filter(int_type != "None")

  drug_hub <- dgi_edges |>
    dplyr::count(drug, sort = TRUE) |>
    dplyr::filter(n >= min_dgidb_targets)

  dgi_edges <- dplyr::semi_join(dgi_edges, drug_hub, by = "drug")

  # Assign DGIdb drug classes by most common interaction type
  drug_int_class <- dgi_edges |>
    dplyr::count(drug, int_type, sort = TRUE) |>
    dplyr::group_by(drug) |>
    dplyr::slice_max(n, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::rename(class = int_type)

  dgi_edges <- dgi_edges |>
    dplyr::left_join(drug_int_class |> dplyr::select(drug, class), by = "drug")

  # ---- Combine edges ----
  curated_edges <- curated_edges |>
    dplyr::select(drug, gene, class, int_type)

  edges_all <- dplyr::bind_rows(
    curated_edges,
    dplyr::anti_join(dgi_edges, curated_edges |> dplyr::select(drug, gene),
                      by = c("drug", "gene"))
  )

  # ---- Interaction style tables ----
  ist <- dgi_interaction_styles

  int_type_colors <- c(
    `Curated MOA` = "#404040",
    stats::setNames(ist$link_color, ist$interaction_type)
  )
  int_type_lty <- c(
    `Curated MOA` = 1L,
    stats::setNames(ist$lty, ist$interaction_type)
  )

  # ---- Adjacency matrix ----
  mat_data <- edges_all |>
    dplyr::count(drug, gene) |>
    tidyr::pivot_wider(names_from = gene, values_from = n, values_fill = 0L) |>
    tibble::column_to_rownames("drug") |>
    as.matrix()

  drug_meta <- edges_all |>
    dplyr::distinct(drug, class) |>
    dplyr::arrange(class, drug)

  drug_order <- drug_meta$drug[drug_meta$drug %in% rownames(mat_data)]
  gene_order <- colnames(mat_data)[order(colSums(mat_data), decreasing = TRUE)]
  mat_ordered <- mat_data[drug_order, gene_order, drop = FALSE]

  # ---- Colours ----
  drug_col_vec <- stats::setNames(
    drug_class_colors[drug_meta$class[match(drug_order, drug_meta$drug)]],
    drug_order
  )
  # Fill missing with grey
  drug_col_vec[is.na(drug_col_vec)] <- "#BBBBBB"

  gene_col_vec <- stats::setNames(rep(gene_color, length(gene_order)), gene_order)
  all_cols <- c(drug_col_vec, gene_col_vec)

  # ---- Link colour & lty matrices ----
  col_mat <- matrix("#F0F0F0", nrow(mat_ordered), ncol(mat_ordered),
                    dimnames = dimnames(mat_ordered))
  lty_mat <- matrix(1L, nrow(mat_ordered), ncol(mat_ordered),
                    dimnames = dimnames(mat_ordered))

  for (i in seq_len(nrow(edges_all))) {
    r <- edges_all$drug[i]
    c_col <- edges_all$gene[i]
    if (r %in% rownames(col_mat) && c_col %in% colnames(col_mat)) {
      it <- edges_all$int_type[i]
      col_mat[r, c_col] <- int_type_colors[it] %||% "#CCCCCC"
      lty_mat[r, c_col] <- int_type_lty[it]    %||% 1L
    }
  }

  # ---- Gaps: group drug classes ----
  class_runs <- drug_meta |>
    dplyr::filter(drug %in% drug_order) |>
    dplyr::mutate(idx = match(drug, drug_order)) |>
    dplyr::group_by(class) |>
    dplyr::summarise(start = min(idx), end = max(idx), .groups = "drop") |>
    dplyr::arrange(start)

  drug_gaps <- rep(1.5, length(drug_order))
  for (i in seq_len(nrow(class_runs))) {
    ei <- class_runs$end[i]
    if (ei < length(drug_order)) drug_gaps[ei] <- 6
  }
  drug_gaps[length(drug_order)] <- 6
  gene_gaps <- rep(1.5, length(gene_order))
  gene_gaps[length(gene_order)] <- 8

  # ---- Draw ----
  .draw_chord <- function() {
    graphics::par(mar = c(1, 1, 1, 1), bg = "white")
    circlize::circos.clear()
    circlize::circos.par(
      gap.after    = c(drug_gaps, gene_gaps),
      start.degree = 85,
      track.margin = c(0.02, 0.08)
    )

    circlize::chordDiagram(
      mat_ordered,
      grid.col            = all_cols,
      col                 = col_mat,
      link.lty            = lty_mat,
      transparency        = transparency,
      directional         = 0,
      annotationTrack     = "grid",
      annotationTrackHeight = 0.03,
      preAllocateTracks = list(
        list(track.height = 0.12),
        list(track.height = 0.12)
      ),
      link.visible = (mat_ordered > 0)
    )

    # Track 1: class colour band
    circlize::circos.track(track.index = 1,
                            panel.fun = function(x, y) {}, bg.border = NA)
    for (i in seq_len(nrow(class_runs))) {
      cls   <- class_runs$class[i]
      s_idx <- class_runs$start[i]
      e_idx <- class_runs$end[i]
      d_list <- drug_order[s_idx:e_idx]
      if (length(d_list) == 0) next
      circlize::highlight.sector(d_list, track.index = 1,
                                  col = drug_class_colors[cls], border = NA)
      mid_drug <- d_list[ceiling(length(d_list) / 2)]
      circlize::set.current.cell(sector.index = mid_drug, track.index = 1)
      circlize::circos.text(
        mean(circlize::get.cell.meta.data("xlim")),
        mean(circlize::get.cell.meta.data("ylim")),
        cls, cex = 0.65, font = 2,
        col = "white",
        facing = "clockwise", niceFacing = TRUE, adj = c(0.5, 0.5)
      )
    }

    # Track 2: sector labels
    circlize::circos.track(
      track.index = 2,
      panel.fun = function(x, y) {
        sn   <- circlize::get.cell.meta.data("sector.index")
        xlim <- circlize::get.cell.meta.data("xlim")
        ylim <- circlize::get.cell.meta.data("ylim")
        is_drug <- sn %in% drug_order
        lbl <- sn
        col <- if (is_drug) "black" else gene_color
        cex <- if (is_drug) 0.9 else 1.1
        fnt <- if (is_drug) 2L  else 1L
        circlize::circos.text(mean(xlim), mean(ylim), lbl,
                               cex = cex, font = fnt, col = col,
                               facing = "clockwise", niceFacing = TRUE,
                               adj = c(0.5, 0.5))
      }, bg.border = NA
    )

    graphics::title(title, cex.main = 2.2, line = 1.5)
    graphics::mtext(paste0(
      length(drug_order), " drugs | ",
      length(gene_order), " biomarkers | ",
      sum(mat_ordered > 0), " edges"
    ), side = 3, line = 0.5, cex = 1.1, col = "grey40")

    # Legend: drug classes
    defined_classes <- intersect(names(drug_class_colors), drug_meta$class)
    if (length(defined_classes) > 0) {
      graphics::legend("bottomleft",
             legend = defined_classes,
             fill   = drug_class_colors[defined_classes],
             border = NA, bty = "n", cex = 0.75,
             title = "Drug class", title.adj = 0, ncol = 2)
    }
  }

  if (!is.null(out_dir)) {
    svg_path <- file.path(out_dir, paste0(base_filename, ".svg"))
    grDevices::svg(svg_path, width = svg_w, height = svg_h, bg = "transparent")
    .draw_chord()
    grDevices::dev.off()
    message("SVG saved: ", svg_path)
  } else {
    .draw_chord()
  }

  invisible(mat_ordered)
}
