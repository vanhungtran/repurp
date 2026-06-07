#' Build an interactive visNetwork drug-gene interaction widget
#'
#' Constructs a self-contained `visNetwork` HTML widget from a standardised
#' DGI data frame. Genes appear as green circles, drugs as yellow squares.
#' Edge colour encodes the interaction type; edge width scales with
#' `interaction_score` when present.
#'
#' @param net_dat A data frame with columns `Gene`, `Drug`, `InteractionType`
#'   (output of [prepare_dgi_data()]).  An optional `interaction_score` column
#'   is used for edge width.
#' @param title Character string shown as the network title. Default `""`.
#' @param gene_color Hex colour for gene nodes. Default `"#009E73"`.
#' @param drug_color Hex colour for drug nodes. Default `"#F0E442"`.
#' @param height Widget height string, e.g. `"800px"`. Default `"800px"`.
#' @param physics_solver visNetwork physics solver. Default `"forceAtlas2Based"`.
#'
#' @return A `visNetwork` htmlwidget object.
#'
#' @export
#' @examples
#' \dontrun{
#' vn <- build_visnetwork(net_dat, title = "DGI — set1")
#' print(vn)
#' }
build_visnetwork <- function(net_dat,
                              title          = "",
                              gene_color     = "#009E73",
                              drug_color     = "#F0E442",
                              height         = "800px",
                              physics_solver = "forceAtlas2Based") {

  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    rlang::abort("Package 'visNetwork' is required. Install with: install.packages('visNetwork')")
  }

  net_dat$InteractionType <- .normalise_interaction(net_dat$InteractionType)
  net_dat$Drug <- ifelse(net_dat$Drug %in% net_dat$Gene,
                          paste0(net_dat$Drug, " [drug]"), net_dat$Drug)

  score_col <- "interaction_score" %in% colnames(net_dat)

  genes_df <- tibble::tibble(
    id    = unique(net_dat$Gene),
    label = unique(net_dat$Gene),
    group = "Gene",
    shape = "dot",
    color = gene_color,
    title = paste0("<b>Gene:</b> ", unique(net_dat$Gene))
  )

  drugs_df <- tibble::tibble(
    id    = unique(net_dat$Drug),
    label = unique(net_dat$Drug),
    group = "Drug",
    shape = "square",
    color = drug_color,
    title = paste0("<b>Drug:</b> ", unique(net_dat$Drug))
  )

  nodes <- dplyr::bind_rows(genes_df, drugs_df)

  edge_col_map <- dgi_edge_colors
  names(edge_col_map) <- toupper(names(edge_col_map))

  edges <- net_dat |>
    dplyr::mutate(
      from  = Gene,
      to    = Drug,
      label = InteractionType,
      color = dplyr::coalesce(
        edge_col_map[toupper(InteractionType)],
        "#999999"
      )
    )

  if (score_col) {
    edges$score_num <- suppressWarnings(as.numeric(net_dat$interaction_score))
    edges$score_num[is.na(edges$score_num)] <- 0
    edges$width <- scales_rescale(edges$score_num, to = c(1, 5))
    edges$title <- paste0(edges$InteractionType,
                           "<br>Score: ", round(edges$score_num, 3))
  } else {
    edges$width <- 1.5
    edges$title <- edges$InteractionType
  }

  edges <- edges[, c("from", "to", "label", "color", "width", "title")]

  vn <- visNetwork::visNetwork(nodes, edges,
                                main   = title,
                                height = height,
                                width  = "100%") |>
    visNetwork::visGroups(groupname = "Gene", color = gene_color,
                           shape = "dot",    size = 20) |>
    visNetwork::visGroups(groupname = "Drug", color = drug_color,
                           shape = "square", size = 15) |>
    visNetwork::visLegend(useGroups = TRUE, position = "right",
                           main = "Node Type") |>
    visNetwork::visEdges(
      smooth = list(enabled = TRUE, type = "continuous"),
      arrows = list(to = list(enabled = TRUE, scaleFactor = 0.5))
    ) |>
    visNetwork::visOptions(
      highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
      nodesIdSelection = TRUE
    ) |>
    visNetwork::visInteraction(
      navigationButtons = TRUE,
      dragNodes  = TRUE,
      dragView   = TRUE,
      zoomView   = TRUE
    ) |>
    visNetwork::visPhysics(
      stabilization = list(iterations = 200),
      solver = physics_solver,
      forceAtlas2Based = list(
        gravitationalConstant = -50,
        centralGravity  = 0.01,
        springLength    = 150,
        springConstant  = 0.08
      )
    ) |>
    visNetwork::visLayout(randomSeed = 42)

  vn
}


#' Save a visNetwork widget to a self-contained HTML file
#'
#' Thin wrapper around [visNetwork::visSave()] that also reports file size.
#'
#' @param vn A `visNetwork` object (from [build_visnetwork()]).
#' @param path Full path to the output HTML file.
#' @param selfcontained Logical; embed all JS/CSS inline? Default `TRUE`.
#'   Set to `FALSE` for large networks to avoid huge files.
#'
#' @return Invisibly returns `path`.
#' @export
#' @examples
#' \dontrun{
#' vn <- build_visnetwork(net_dat)
#' save_visnetwork(vn, "results/Fig15A.html")
#' }
save_visnetwork <- function(vn, path, selfcontained = TRUE) {
  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    rlang::abort("Package 'visNetwork' is required.")
  }
  visNetwork::visSave(vn, path, selfcontained = selfcontained)
  message("Saved: ", path,
          " (", round(file.info(path)$size / 1024, 1), " KB)")
  invisible(path)
}

# thin shim so we don't hard-depend on {scales}
scales_rescale <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
  if (diff(from) == 0) return(rep(mean(to), length(x)))
  (x - from[1]) / diff(from) * diff(to) + to[1]
}
