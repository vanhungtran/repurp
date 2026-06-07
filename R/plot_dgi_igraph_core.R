#' Plot a core drug–gene interaction network with igraph (high-connectivity genes only)
#'
#' Builds on [plot_dgi_igraph()] but filters to genes with at least
#' `min_gene_edges` drug connections before plotting, reducing visual clutter
#' for large networks. Uses the same igraph + qgraph FR layout engine.
#'
#' @inheritParams plot_dgi_igraph
#' @param min_gene_edges Minimum number of drug connections for a gene to be
#'   included. Default `3`.
#'
#' @return Invisibly returns a named list with elements `graph` (the `igraph`
#'   object) and `layout` (the numeric layout matrix).
#' @export
#' @examples
#' \dontrun{
#' net_dat <- prepare_dgi_data("interactions.tsv", my_genes)
#' plot_dgi_igraph_core(net_dat, title = "Core DGI Network", out_dir = "results/")
#' }
plot_dgi_igraph_core <- function(net_dat,
                                  title          = "Core Drug–Gene Interaction Network",
                                  out_dir        = NULL,
                                  base_filename  = paste0("DGI_core_", Sys.Date()),
                                  gene_color     = "#009E73",
                                  drug_color     = "#F0E442",
                                  vertex_size    = 3,
                                  label_cex      = 0.5,
                                  area_factor    = 2.05,
                                  repulse_factor = 3,
                                  min_gene_edges = 3L,
                                  svg_w = 22, svg_h = 20,
                                  png_w = 2200, png_h = 2000,
                                  save_svg = !is.null(out_dir),
                                  save_png = !is.null(out_dir)) {

  net_dat <- .prep_net_dat(net_dat)

  # Filter to core genes
  gene_degree <- net_dat |>
    dplyr::count(Gene, name = "n_drugs")

  core_genes <- gene_degree |>
    dplyr::filter(n_drugs >= min_gene_edges) |>
    dplyr::pull(Gene)

  if (length(core_genes) == 0) {
    rlang::abort(sprintf(
      "No genes with >= %d drug connections. Try lowering `min_gene_edges`.",
      min_gene_edges
    ))
  }

  net_dat <- net_dat |> dplyr::filter(Gene %in% core_genes)

  message(sprintf(
    "Core network: %d genes, %d drugs, %d edges",
    length(unique(net_dat$Gene)),
    length(unique(net_dat$Drug)),
    nrow(net_dat)
  ))

  edges <- net_dat
  colnames(edges)[colnames(edges) == "Gene"] <- "from"
  colnames(edges)[colnames(edges) == "Drug"] <- "to"

  # Disambiguate drug names that collide with gene names
  collide <- edges$to %in% edges$from
  edges$to[collide] <- paste0(edges$to[collide], "_drug")

  nodes <- data.frame(
    id    = c(unique(edges$from), unique(edges$to)),
    label = c(unique(edges$from), unique(edges$to)),
    type  = c(rep("gene", length(unique(edges$from))),
              rep("drug", length(unique(edges$to)))),
    stringsAsFactors = FALSE
  )

  net <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

  coul  <- grDevices::rainbow(nlevels(as.factor(edges$InteractionType)))
  my_ec <- coul[as.numeric(as.factor(edges$InteractionType))]
  coul1 <- c(gene = gene_color, drug = drug_color)
  my_nc <- coul1[nodes$type]

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
         vertex.color       = my_nc,
         vertex.label.cex   = label_cex,
         vertex.frame.color = "transparent",
         vertex.label.color  = "black",
         vertex.label.family = "sans",
         vertex.label.font   = 1L,
         vertex.label.dist   = 0,
         vertex.label.degree = 0,
         edge.color    = my_ec,
         edge.width    = 1,
         edge.arrow.size  = 1,
         edge.arrow.width = 1,
         edge.lty     = "solid",
         edge.curved  = 0.3)

    graphics::text(-1, 1.2, title, col = "black", cex = 1.5)

    graphics::legend(x = 1.1, y = 0.2,
           legend = paste(levels(as.factor(edges$InteractionType)), " "),
           col = coul, bty = "n", pch = 20, pt.cex = 2, cex = 1,
           text.col = "black", horiz = FALSE,
           title = "Interaction type")

    graphics::legend(x = 1.1, y = 0.4,
           legend = paste(levels(as.factor(nodes$type)), " "),
           col = coul1, bty = "n", pch = 20, pt.cex = 2, cex = 1,
           text.col = "black", horiz = FALSE,
           title = "Node type")
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
