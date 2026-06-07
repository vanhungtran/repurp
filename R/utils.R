#' Edge colour palette keyed by DGIdb interaction type
#'
#' A named character vector mapping DGIdb interaction type labels to hex colours.
#' Used consistently across all Repurp plot functions.
#'
#' @export
#' @examples
#' dgi_edge_colors["INHIBITOR"]
dgi_edge_colors <- c(
  INHIBITOR      = "#E41A1C",
  AGONIST        = "#377EB8",
  ANTAGONIST     = "#FF7F00",
  ACTIVATOR      = "#4DAF4A",
  BINDER         = "#984EA3",
  BLOCKER        = "#A65628",
  MODULATOR      = "#F781BF",
  SUBSTRATE      = "#66C2A5",
  `Inverse Agonist` = "#FB8072",
  Antibody       = "#8DD3C7",
  Cleavage       = "#BEBADA",
  Other          = "#BDBDBD",
  None           = "#999999",
  `Curated MOA`  = "#404040"
)

#' Node colour palette keyed by node type
#'
#' @export
#' @examples
#' dgi_node_colors["Gene"]
dgi_node_colors <- c(
  Gene    = "#009E73",
  Drug    = "#F0E442",
  Pathway = "#6A3D9A"
)

#' Full interaction-style table for chord and sector colouring
#'
#' A tibble with columns `interaction_type`, `sector_color`, `link_color`, and `lty`.
#'
#' @export
dgi_interaction_styles <- tibble::tribble(
  ~interaction_type,  ~sector_color, ~link_color, ~lty,
  "Inhibitor",        "#FB9A99",    "#FBB4AE",   2L,
  "Agonist",          "#A6CEE3",    "#B3CDE3",   3L,
  "Antagonist",       "#FDBF6F",    "#FED9A6",   4L,
  "Binder",           "#CAB2D6",    "#DECBE4",   5L,
  "Blocker",          "#B15928",    "#CCEBC5",   6L,
  "Activator",        "#B2DF8A",    "#E5F5E0",   2L,
  "Modulator",        "#FCCDE5",    "#FFF7BC",   3L,
  "Substrate",        "#80CDC1",    "#C2E699",   4L,
  "Antibody",         "#8DD3C7",    "#A6D854",   5L,
  "Cleavage",         "#BEBADA",    "#D4B9DA",   6L,
  "Inverse Agonist",  "#FB8072",    "#FDD0A2",   2L,
  "Other",            "#BDBDBD",    "#D9D9D9",   1L
)

# Internal helper: normalise InteractionType strings
.normalise_interaction <- function(x) {
  dplyr::case_when(
    is.na(x) | x == "NULL" | x == "" ~ "None",
    TRUE ~ x
  )
}

# Internal: check ggraph/tidygraph/ggplot2 availability
.check_ggraph <- function() {
  for (pkg in c("ggraph", "tidygraph", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      rlang::abort(sprintf(
        "Package '%s' is required. Install with: install.packages('%s')", pkg, pkg
      ))
    }
  }
}

# Internal: check optional package availability
.check_optional_pkg <- function(pkg, caller = "") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    rlang::abort(sprintf(
      "Package '%s' is required for %s. Install with: install.packages('%s')",
      pkg, caller, pkg
    ))
  }
}

# Internal: save ggplot to PNG + SVG
.save_ggplot <- function(p, out_dir, base_filename, width, height, dpi) {
  if (is.null(out_dir)) {
    print(p)
    return(invisible())
  }
  png_path <- file.path(out_dir, paste0(base_filename, ".png"))
  ggplot2::ggsave(png_path, p, width = width, height = height, dpi = dpi, bg = "white")
  message("PNG saved: ", png_path)

  svg_path <- file.path(out_dir, paste0(base_filename, ".svg"))
  ggplot2::ggsave(svg_path, p, width = width, height = height, bg = "white")
  message("SVG saved: ", svg_path)
}
