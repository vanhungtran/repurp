#' Plot a drug-gene interaction network with igraph + Fruchterman-Reingold layout
#'
#' Draws a bipartite drug–gene network using an `igraph` graph object laid out with
#' `qgraph::qgraph.layout.fruchtermanreingold()`. Optionally saves the result as
#' SVG (vector, poster-quality) and/or PNG (raster preview).
#'
#' @param net_dat A data frame with columns `Gene`, `Drug`, `InteractionType` (and
#'   optionally `Count` / other columns). Typically the output of [prepare_dgi_data()].
#' @param title Character string for the network title (rendered via `text()`).
#' @param out_dir Directory path for saved files. If `NULL` (default) only the plot is
#'   drawn to the current device; no files are written.
#' @param base_filename Base name (without extension) for saved files. Defaults to
#'   `paste0("DGI_network_", Sys.Date())`.
#' @param gene_color Hex colour for gene nodes. Default `"#009E73"`.
#' @param drug_color Hex colour for drug nodes. Default `"#F0E442"`.
#' @param vertex_size Numeric size passed to `igraph::plot.igraph`. Default `4`.
#' @param label_cex Character expansion for node labels. Default `0.5`.
#' @param area_factor Exponent for the FR layout `area` parameter (`vcount^area_factor`).
#'   Default `2.05`.
#' @param repulse_factor Exponent for the FR `repulse.rad` parameter. Default `3`.
#' @param svg_w,svg_h Width/height in inches for SVG output. Defaults `20` / `18`.
#' @param png_w,png_h Width/height in pixels for PNG output. Defaults `2000` / `1750`.
#' @param save_svg Logical; write an SVG file? Default `TRUE` when `out_dir` is given.
#' @param save_png Logical; write a PNG file? Default `TRUE` when `out_dir` is given.
#'
#' @return Invisibly returns a named list with elements `graph` (the `igraph` object)
#'   and `layout` (the numeric layout matrix).
#'
#' @export
#' @examples
#' \dontrun{
#' net_dat <- prepare_dgi_data("interactions.tsv", my_genes)
#' plot_dgi_igraph(net_dat, title = "My DGI Network", out_dir = "results/")
#' }
plot_dgi_igraph <- function(net_dat,
                             title          = "Drug-Gene Interaction Network",
                             out_dir        = NULL,
                             base_filename  = paste0("DGI_network_", Sys.Date()),
                             gene_color     = "#009E73",
                             drug_color     = "#F0E442",
                             vertex_size    = 4,
                             label_cex      = 0.5,
                             area_factor    = 2.05,
                             repulse_factor = 3,
                             svg_w = 20, svg_h = 18,
                             png_w = 2000, png_h = 1750,
                             save_svg = !is.null(out_dir),
                             save_png = !is.null(out_dir)) {

  net_dat <- .prep_net_dat(net_dat)

  edges <- net_dat
  colnames(edges)[colnames(edges) == "Gene"] <- "from"
  colnames(edges)[colnames(edges) == "Drug"] <- "to"

  nodes <- data.frame(
    id    = c(unique(net_dat$Gene), unique(net_dat$Drug)),
    label = c(unique(net_dat$Gene), unique(net_dat$Drug)),
    type  = c(rep("gene", length(unique(net_dat$Gene))),
              rep("drug", length(unique(net_dat$Drug)))),
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

  # Always draw to current device when called interactively
  if (is.null(out_dir) || (!save_svg && !save_png)) {
    .do_plot()
  }

  message("  Genes: ", length(unique(net_dat$Gene)),
          "  Drugs: ", length(unique(net_dat$Drug)),
          "  Edges: ", nrow(edges))

  invisible(list(graph = net, layout = l))
}

# ---- internal ----
.prep_net_dat <- function(net_dat) {
  net_dat$InteractionType <- .normalise_interaction(net_dat$InteractionType)
  net_dat$Count <- as.integer(table(net_dat$Gene)[net_dat$Gene])
  net_dat
}
