#' Build an interactive tripartite visNetwork (Gene + Drug + Pathway)
#'
#' Constructs a self-contained `visNetwork` HTML widget with three node types:
#' Gene (green circle), Drug (yellow square), and Pathway (purple diamond).
#' Pathway edges typically come from STRING enrichment (KEGG, WikiPathways,
#' GO Process). Node size scales with connectivity; edge colour encodes
#' interaction type or pathway category.
#'
#' @inheritParams build_visnetwork
#' @param pathway_edges A data frame with columns `from` (gene), `to` (pathway
#'   description), `edge_type` (category e.g. "KEGG", "WikiPathways", "GO_Process"),
#'   and optionally `fdr` (numeric). Pass `NULL` to omit pathways.
#' @param pathway_color Hex colour for pathway nodes. Default `"#6A3D9A"`.
#' @param pathway_shape visNetwork shape for pathway nodes. Default `"diamond"`.
#' @param height Widget height string. Default `"900px"`.
#'
#' @return A `visNetwork` htmlwidget object.
#' @export
#' @examples
#' \dontrun{
#' # requires rbioapi for STRING enrichment
#' enrich <- rbioapi::rba_string_enrichment(ids = my_genes, species = 9606)
#' pw_edges <- dplyr::bind_rows(
#'   enrich[["KEGG"]]         |> dplyr::mutate(category = "KEGG"),
#'   enrich[["WikiPathways"]] |> dplyr::mutate(category = "WikiPathways"),
#'   enrich[["Process"]]      |> dplyr::mutate(category = "GO_Process")
#' ) |> dplyr::filter(fdr < 0.05) |>
#'   tidyr::unnest(inputGenes) |>
#'   dplyr::transmute(from = inputGenes, to = description, edge_type = category, fdr = fdr)
#'
#' vn <- build_visnetwork_tripartite(net_dat, pw_edges, title = "Tripartite Network")
#' print(vn)
#' }
build_visnetwork_tripartite <- function(net_dat,
                                         pathway_edges   = NULL,
                                         title           = "Drug–Gene–Pathway Tripartite Network",
                                         gene_color      = "#009E73",
                                         drug_color      = "#F0E442",
                                         pathway_color   = "#6A3D9A",
                                         pathway_shape   = "diamond",
                                         height          = "900px",
                                         physics_solver  = "forceAtlas2Based") {

  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    rlang::abort("Package 'visNetwork' is required. Install with: install.packages('visNetwork')")
  }

  net_dat$InteractionType <- .normalise_interaction(net_dat$InteractionType)
  net_dat$Drug <- ifelse(net_dat$Drug %in% net_dat$Gene,
                          paste0(net_dat$Drug, " [drug]"), net_dat$Drug)

  score_col <- "interaction_score" %in% colnames(net_dat)

  # ---- Build DGI edges ----
  edge_col_map <- dgi_edge_colors
  names(edge_col_map) <- toupper(names(edge_col_map))

  dgi_edges <- net_dat |>
    dplyr::mutate(
      from  = Gene,
      to    = Drug,
      label = InteractionType,
      color = dplyr::coalesce(
        edge_col_map[toupper(InteractionType)],
        "#999999"
      ),
      edge_category = "drug_gene"
    )

  if (score_col) {
    dgi_edges$score_num <- suppressWarnings(as.numeric(net_dat$interaction_score))
    dgi_edges$score_num[is.na(dgi_edges$score_num)] <- 0
    dgi_edges$width <- scales_rescale(dgi_edges$score_num, to = c(1, 5))
    dgi_edges$title <- paste0(dgi_edges$InteractionType,
                               "<br>Score: ", round(dgi_edges$score_num, 3))
  } else {
    dgi_edges$width <- 1.5
    dgi_edges$title <- dgi_edges$InteractionType
  }

  # ---- Build pathway edges ----
  path_ids <- character(0)
  if (!is.null(pathway_edges)) {
    path_ids <- unique(pathway_edges$to)

    gp_edges <- pathway_edges |>
      dplyr::mutate(
        from  = from,
        to    = to,
        label = edge_type,
        color = dplyr::case_when(
          edge_type == "KEGG"          ~ "#1F78B4",
          edge_type == "WikiPathways"  ~ "#33A02C",
          edge_type == "GO_Process"    ~ "#FB9A99",
          edge_type == "GO Process"    ~ "#FB9A99",
          TRUE                         ~ "#6A3D9A"
        ),
        width = if ("fdr" %in% colnames(pathway_edges)) {
          scales_rescale(1 / (as.numeric(fdr) + 1e-10), to = c(1, 4))
        } else {
          2
        },
        title = paste0(to, "<br>Category: ", edge_type),
        edge_category = "gene_pathway"
      ) |>
      dplyr::select(from, to, label, color, width, title, edge_category)

    all_edges <- dplyr::bind_rows(
      dgi_edges |> dplyr::select(from, to, label, color, width, title, edge_category),
      gp_edges
    )
  } else {
    all_edges <- dgi_edges |>
      dplyr::select(from, to, label, color, width, title, edge_category)
  }

  # ---- Build nodes ----
  genes_df <- tibble::tibble(
    id    = unique(net_dat$Gene),
    label = unique(net_dat$Gene),
    group = "Gene",
    shape = "dot",
    color = gene_color,
    size  = 20,
    title = paste0("<b>Gene:</b> ", unique(net_dat$Gene))
  )

  drugs_df <- tibble::tibble(
    id    = unique(net_dat$Drug),
    label = unique(net_dat$Drug),
    group = "Drug",
    shape = "square",
    color = drug_color,
    size  = 15,
    title = paste0("<b>Drug:</b> ", unique(net_dat$Drug))
  )

  nodes <- dplyr::bind_rows(genes_df, drugs_df)

  if (!is.null(pathway_edges)) {
    pw_summary <- pathway_edges |>
      dplyr::group_by(to, edge_type) |>
      dplyr::summarise(
        n_genes = dplyr::n(),
        min_fdr = if ("fdr" %in% colnames(pathway_edges))
          min(as.numeric(fdr), na.rm = TRUE) else NA_real_,
        .groups = "drop"
      )

    pathway_nodes <- tibble::tibble(
      id    = pw_summary$to,
      label = pw_summary$to,
      group = "Pathway",
      shape = pathway_shape,
      color = dplyr::case_when(
        pw_summary$edge_type == "KEGG"          ~ "#1F78B4",
        pw_summary$edge_type == "WikiPathways"  ~ "#33A02C",
        pw_summary$edge_type %in% c("GO_Process", "GO Process") ~ "#FB9A99",
        TRUE                                    ~ pathway_color
      ),
      size  = 20 + scales_rescale(pw_summary$n_genes, to = c(0, 15)),
      title = paste0("<b>Pathway:</b> ", pw_summary$to,
                     "<br>Category: ", pw_summary$edge_type,
                     "<br>Genes: ", pw_summary$n_genes,
                     if (!is.na(pw_summary$min_fdr[1]))
                       paste0("<br>FDR: ", formatC(pw_summary$min_fdr,
                                                    format = "e", digits = 2))
                     else "")
    )

    nodes <- dplyr::bind_rows(nodes, pathway_nodes) |>
      dplyr::distinct(id, .keep_all = TRUE)
  }

  # ---- Build visNetwork ----
  n_genes  <- nrow(genes_df)
  n_drugs  <- nrow(drugs_df)
  n_path   <- length(path_ids)
  n_edges  <- nrow(all_edges)

  vn <- visNetwork::visNetwork(
    nodes, all_edges,
    main    = title,
    submain = paste0(
      "Genes: ", n_genes, " | Drugs: ", n_drugs,
      if (n_path > 0) paste0(" | Pathways: ", n_path) else "",
      " | Edges: ", n_edges
    ),
    height  = height,
    width   = "100%"
  ) |>
    visNetwork::visGroups(groupname = "Gene", color = gene_color,
                           shape = "dot",     size = 20) |>
    visNetwork::visGroups(groupname = "Drug", color = drug_color,
                           shape = "square",  size = 15)

  if (n_path > 0) {
    vn <- vn |>
      visNetwork::visGroups(groupname = "Pathway", color = pathway_color,
                             shape = pathway_shape, size = 25)
  }

  vn <- vn |>
    visNetwork::visLegend(useGroups = TRUE, position = "right",
                           main = "Node Type", ncol = 1) |>
    visNetwork::visEdges(
      smooth = list(enabled = TRUE, type = "continuous"),
      arrows = list(to = list(enabled = TRUE, scaleFactor = 0.5))
    ) |>
    visNetwork::visOptions(
      highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
      nodesIdSelection = list(enabled = TRUE,
                               values = sort(nodes$id),
                               style = "width: 300px;"),
      selectedBy = list(variable = "group", multiple = TRUE,
                         main = "Filter by node type")
    ) |>
    visNetwork::visInteraction(
      navigationButtons = TRUE,
      dragNodes  = TRUE,
      dragView   = TRUE,
      zoomView   = TRUE
    ) |>
    visNetwork::visPhysics(
      stabilization = list(iterations = 300),
      solver = physics_solver,
      forceAtlas2Based = list(
        gravitationalConstant = -80,
        centralGravity  = 0.005,
        springLength    = 200,
        springConstant  = 0.05
      )
    ) |>
    visNetwork::visLayout(randomSeed = 42)

  vn
}
