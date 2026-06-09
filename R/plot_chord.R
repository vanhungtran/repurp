#' Draw a circlize chord diagram of drug–gene interactions
#'
#' Renders a circular chord diagram (via the `circlize` package) showing
#' drug–biomarker interaction edges. Drug sectors are coloured by drug class and
#' grouped into a continuous colour band; gene sectors are coloured green.
#' Link ribbon colour encodes interaction type.
#'
#' @param edges A data frame with columns `drug`, `gene`, `class` (drug class
#'   string), and `int_type` (interaction-type string).  Typically built by
#'   combining [repurp_ad_edges()] with DGIdb edges; see the package vignette.
#' @param drug_class_colors Named character vector mapping drug-class strings to
#'   hex colours.  Defaults to [repurp_class_colors()].
#' @param gene_color Hex colour for gene sectors. Default `"#009E73"`.
#' @param title Plot title. Default `"Drug–Gene Chord Diagram"`.
#' @param out_dir Directory for saved SVG file. `NULL` = draw to current device only.
#' @param base_filename Base filename without extension.
#'   Default `paste0("DGI_chord_", Sys.Date())`.
#' @param svg_w,svg_h SVG dimensions in inches. Defaults `32` / `32`.
#' @param transparency Chord transparency (0–1). Default `0.15`.
#' @param min_drug_targets Integer; drugs with fewer gene targets than this are
#'   dropped to reduce clutter. Default `1L` (keep all).
#'
#' @return Invisibly returns the adjacency matrix used for plotting.
#' @export
#' @examples
#' \dontrun{
#' edges <- dplyr::bind_rows(repurp_ad_edges(), my_dgidb_edges)
#' plot_chord_dgi(edges, out_dir = "results/")
#' }
plot_chord_dgi <- function(edges,
                            drug_class_colors = repurp_class_colors(),
                            gene_color        = "#009E73",
                            title             = "Drug–Gene Chord Diagram",
                            out_dir           = NULL,
                            base_filename     = paste0("DGI_chord_", Sys.Date()),
                            svg_w             = 32,
                            svg_h             = 32,
                            transparency      = 0.15,
                            min_drug_targets  = 1L) {

  if (!requireNamespace("circlize", quietly = TRUE)) {
    rlang::abort("Package 'circlize' is required. Install with: install.packages('circlize')")
  }

  # Drop low-connectivity drugs
  if (min_drug_targets > 1L) {
    hub <- edges |>
      dplyr::count(drug) |>
      dplyr::filter(n >= min_drug_targets)
    edges <- dplyr::semi_join(edges, hub, by = "drug")
  }

  # ---- Adjacency matrix ----
  mat_data <- edges |>
    dplyr::count(drug, gene) |>
    tidyr::pivot_wider(names_from = gene, values_from = n, values_fill = 0L) |>
    tibble::column_to_rownames("drug") |>
    as.matrix()

  # Order: sort drugs by class, then name; genes by degree descending
  drug_meta <- edges |>
    dplyr::distinct(drug, class) |>
    dplyr::arrange(class, drug)

  drug_order <- drug_meta$drug[drug_meta$drug %in% rownames(mat_data)]
  gene_order  <- colnames(mat_data)[order(colSums(mat_data), decreasing = TRUE)]
  mat_ordered <- mat_data[drug_order, gene_order, drop = FALSE]

  # ---- Sector colours ----
  drug_col_vec <- stats::setNames(
    drug_class_colors[drug_meta$class[match(drug_order, drug_meta$drug)]],
    drug_order
  )
  gene_col_vec <- stats::setNames(rep(gene_color, length(gene_order)), gene_order)
  all_cols <- c(drug_col_vec, gene_col_vec)

  # ---- Interaction-type link colours ----
  ist <- dgi_interaction_styles
  int_type_colors <- c(
    `Curated MOA` = "#404040",
    stats::setNames(ist$link_color, ist$interaction_type)
  )
  int_type_lty <- c(
    `Curated MOA` = 1L,
    stats::setNames(ist$lty, ist$interaction_type)
  )

  col_mat <- matrix("#F0F0F0", nrow(mat_ordered), ncol(mat_ordered),
                    dimnames = dimnames(mat_ordered))
  lty_mat <- matrix(1L, nrow(mat_ordered), ncol(mat_ordered),
                    dimnames = dimnames(mat_ordered))

  for (i in seq_len(nrow(edges))) {
    r <- edges$drug[i];  c_col <- edges$gene[i]
    if (r %in% rownames(col_mat) && c_col %in% colnames(col_mat)) {
      it <- edges$int_type[i]
      col_mat[r, c_col] <- int_type_colors[it] %||% "#CCCCCC"
      lty_mat[r, c_col] <- int_type_lty[it]    %||% 1L
    }
  }

  # ---- Class runs for gap and band ----
  class_runs <- drug_meta |>
    dplyr::filter(drug %in% drug_order) |>
    dplyr::mutate(idx = match(drug, drug_order)) |>
    dplyr::group_by(class) |>
    dplyr::summarise(start = min(idx), end = max(idx), .groups = "drop") |>
    dplyr::arrange(start)

  drug_gaps <- rep(1.5, length(drug_order))
  for (i in seq_len(nrow(class_runs))) {
    ei <- class_runs$end[i]
    if (ei < length(drug_order)) drug_gaps[ei] <- 4
  }
  drug_gaps[length(drug_order)] <- 5
  gene_gaps <- rep(1.5, length(gene_order))
  gene_gaps[length(gene_order)] <- 6

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
      cls     <- class_runs$class[i]
      s_idx   <- class_runs$start[i]
      e_idx   <- class_runs$end[i]
      d_list  <- drug_order[s_idx:e_idx]
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
        lbl  <- sn
        col  <- if (is_drug) "black" else gene_color
        cex  <- if (is_drug) 0.9 else 1.1
        fnt  <- if (is_drug) 2L  else 1L
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
    graphics::legend("bottomleft",
           legend = defined_classes,
           fill   = drug_class_colors[defined_classes],
           border = NA, bty = "n", cex = 0.75,
           title = "Drug class", title.adj = 0, ncol = 2)

    circlize::circos.clear()
  }

  if (!is.null(out_dir)) {
    svg_path <- file.path(out_dir, paste0(base_filename, ".svg"))
    grDevices::svg(svg_path, width = svg_w, height = svg_h)
    .draw_chord()
    grDevices::dev.off()
    message("SVG saved: ", svg_path)
  } else {
    .draw_chord()
  }

  invisible(mat_ordered)
}


#' Default drug-class colour palette for chord diagrams
#'
#' @return A named character vector.
#' @export
repurp_class_colors <- function() {
  c(
    "Corticosteroid"          = "#E41A1C",
    "Calcineurin inhibitor"   = "#FF7F00",
    "PDE4 inhibitor"          = "#377EB8",
    "JAK inhibitor"           = "#4DAF4A",
    "Biologic"                = "#984EA3",
    "AhR agonist"             = "#F781BF",
    "Immunosuppressant"       = "#A65628",
    "Inhibitor"               = "#FB9A99",
    "Agonist"                 = "#A6CEE3",
    "Antagonist"              = "#FDBF6F",
    "Binder"                  = "#CAB2D6",
    "Activator"               = "#B2DF8A",
    "Other"                   = "#BDBDBD"
  )
}

# null-coalesce helper (not importing rlang::`%||%` to keep deps light)
`%||%` <- function(a, b) if (!is.null(a) && !is.na(a)) a else b
