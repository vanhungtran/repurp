#' Plot a tripartite Drug–Gene–Pathway network with igraph + FR layout
#'
#' Builds a three-node-type network (Gene, Drug, Pathway) using `igraph` and
#' `qgraph` Fruchterman-Reingold layout. Pathway edges (e.g. from STRING
#' enrichment via \code{rbioapi::rba_string_enrichment()}) are passed as
#' a separate data frame.
#'
#' @param dgi_edges A data frame with columns `from` (gene), `to` (drug),
#'   `edge_type` (interaction type string), and optionally `source` (e.g. "DGIdb").
#' @param pathway_edges A data frame with columns `from` (gene), `to` (pathway
#'   description), `edge_type` (category string, e.g. "KEGG"), and optionally
#'   `source` (e.g. "STRING"). Pass `NULL` (default) to omit pathways.
#' @param title Plot title.
#' @param out_dir Directory for saved files. `NULL` = draw to current device only.
#' @param base_filename Base filename without extension.
#'   Default `paste0("DGI_tripartite_", Sys.Date())`.
#' @param gene_color,drug_color,pathway_color Hex colours for the three node types.
#'   Defaults `"#009E73"`, `"#F0E442"`, `"#6A3D9A"`.
#' @param vertex_size Numeric vertex size. Default `3`.
#' @param label_cex Character expansion for labels. Default `0.45`.
#' @param area_factor,repulse_factor FR layout exponents. Defaults `2.05` / `3`.
#' @param svg_w,svg_h,png_w,png_h SVG/PNG dimensions.
#' @param save_svg,save_png Logical; write output files?
#'
#' @return Invisibly returns a named list with elements `graph` (the `igraph`
#'   object) and `layout` (the numeric layout matrix).
#' @export
#' @examples
#' \dontrun{
#' # requires rbioapi for STRING enrichment
#' enrich <- rbioapi::rba_string_enrichment(ids = my_genes, species = 9606)
#' pw_edges <- dplyr::bind_rows(
#'   enrich[["KEGG"]]         |> dplyr::mutate(category = "KEGG"),
#'   enrich[["WikiPathways"]] |> dplyr::mutate(category = "WikiPathways"),
#'   enrich[["Process"]]      |> dplyr::mutate(category = "GO Process")
#' ) |> dplyr::filter(fdr < 0.05) |>
#'   tidyr::unnest(inputGenes) |>
#'   dplyr::transmute(from = inputGenes, to = description, edge_type = category)
#'
#' plot_dgi_tripartite_igraph(
#'   dgi_edges    = dplyr::transmute(net_dat, from = Gene, to = Drug, edge_type = InteractionType),
#'   pathway_edges = pw_edges,
#'   out_dir = "results/"
#' )
#' }
plot_dgi_tripartite_igraph <- function(dgi_edges,
                                        pathway_edges   = NULL,
                                        title           = "Drug–Gene–Pathway Network",
                                        out_dir         = NULL,
                                        base_filename   = paste0("DGI_tripartite_", Sys.Date()),
                                        gene_color      = "#009E73",
                                        drug_color      = "#F0E442",
                                        pathway_color   = "#6A3D9A",
                                        vertex_size     = 3,
                                        label_cex       = 0.45,
                                        area_factor     = 2.05,
                                        repulse_factor  = 3,
                                        svg_w = 24, svg_h = 22,
                                        png_w = 2400, png_h = 2200,
                                        save_svg = !is.null(out_dir),
                                        save_png = !is.null(out_dir)) {

  dgi_edges$edge_type <- .normalise_interaction(dgi_edges$edge_type)

  gene_ids <- unique(dgi_edges$from)
  drug_ids <- unique(dgi_edges$to)

  all_edges <- dgi_edges |>
    dplyr::mutate(source = if ("source" %in% colnames(dgi_edges)) source else "DGIdb")

  path_ids <- character(0)
  if (!is.null(pathway_edges)) {
    path_ids <- unique(pathway_edges$to)
    gp <- pathway_edges |>
      dplyr::mutate(source = if ("source" %in% colnames(pathway_edges)) source else "STRING")
    all_edges <- dplyr::bind_rows(all_edges, gp)
  }

  # ---- Build nodes ----
  nodes <- data.frame(
    id    = c(gene_ids, drug_ids, path_ids),
    label = c(gene_ids, drug_ids, path_ids),
    type  = c(rep("Gene", length(gene_ids)),
              rep("Drug", length(drug_ids)),
              rep("Pathway", length(path_ids))),
    stringsAsFactors = FALSE
  )

  # Disambiguate drug names that collide with gene names
  collide <- all_edges$to %in% all_edges$from
  all_edges$to[collide] <- paste0(all_edges$to[collide], "_drug")

  net <- igraph::graph_from_data_frame(all_edges, vertices = nodes, directed = FALSE)

  # ---- Colour palettes ----
  edge_palette <- c(
    INHIBITOR  = "#E41A1C", AGONIST    = "#377EB8",
    ANTAGONIST = "#FF7F00", BINDER     = "#984EA3",
    BLOCKER    = "#A65628", ACTIVATOR  = "#4DAF4A",
    MODULATOR  = "#F781BF", SUBSTRATE  = "#66C2A5",
    None       = "#CCCCCC",
    KEGG       = "#1F78B4", `WikiPathways` = "#33A02C",
    `GO Process` = "#FB9A99"
  )

  edge_types_present <- unique(all_edges$edge_type)
  edge_palette <- edge_palette[names(edge_palette) %in% edge_types_present]

  my_color <- edge_palette[all_edges$edge_type]
  my_color[is.na(my_color)] <- "#CCCCCC"

  node_palette <- c(Gene = gene_color, Drug = drug_color, Pathway = pathway_color)
  my_color1 <- node_palette[nodes$type]

  # ---- qgraph FR layout ----
  e <- igraph::as_edgelist(net, names = FALSE)
  l <- qgraph::qgraph.layout.fruchtermanreingold(
    e, vcount = igraph::vcount(net),
    area       = 2 * (igraph::vcount(net)^area_factor),
    repulse.rad = (igraph::vcount(net)^repulse_factor)
  )

  .do_plot <- function() {
    graphics::par(mar = c(0.5, 0.5, 0.5, 0.5))
    plot(net, layout = l,
         vertex.size        = vertex_size,
         vertex.color       = my_color1,
         vertex.label.cex   = label_cex,
         vertex.frame.color = "transparent",
         vertex.label.color  = "black",
         vertex.label.family = "sans",
         vertex.label.font   = 1L,
         vertex.label.dist   = 0,
         vertex.label.degree = 0,
         edge.color    = my_color,
         edge.width    = 1,
         edge.arrow.size  = 0.8,
         edge.arrow.width = 0.8,
         edge.lty     = "solid",
         edge.curved  = 0.3)

    graphics::text(-1, 1.2, title, col = "black", cex = 1.5)
    graphics::mtext(paste0(
      length(gene_ids), " genes | ", length(drug_ids), " drugs | ",
      if (!is.null(pathway_edges)) paste0(length(path_ids), " pathways | ") else "",
      igraph::ecount(net), " edges"
    ), side = 3, line = -1.5, cex = 0.7, col = "grey40")

    # Edge type legend
    if (length(edge_palette) > 0) {
      graphics::legend(x = 1.1, y = 0.5,
             legend = paste(names(edge_palette), " "),
             col = edge_palette, bty = "n", pch = 20, pt.cex = 1.5, cex = 0.7,
             text.col = "black", horiz = FALSE, ncol = 2)
    }

    # Node type legend
    graphics::legend(x = 1.1, y = 0.25,
           legend = paste(names(node_palette), " "),
           col = node_palette, bty = "n", pch = 20, pt.cex = 1.8, cex = 0.85,
           text.col = "black", horiz = FALSE)
  }

  if (save_svg && !is.null(out_dir)) {
    svg_path <- file.path(out_dir, paste0(base_filename, ".svg"))
    grDevices::svg(svg_path, width = svg_w, height = svg_h, bg = "transparent")
    .do_plot()
    grDevices::dev.off()
    message("SVG saved: ", svg_path)
  }

  if (save_png && !is.null(out_dir)) {
    png_path <- file.path(out_dir, paste0(base_filename, ".png"))
    grDevices::png(png_path, width = png_w, height = png_h, res = 96,
                   bg = "transparent")
    .do_plot()
    grDevices::dev.off()
    message("PNG saved: ", png_path)
  }

  if (is.null(out_dir) || (!save_svg && !save_png)) {
    .do_plot()
  }

  invisible(list(graph = net, layout = l))
}
