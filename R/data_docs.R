#' Curated reference table of 18 established atopic-dermatitis drugs
#'
#' A tibble with one row per drug, covering all currently approved or
#' commonly used systemic and topical AD treatments.
#'
#' @format A tibble with 18 rows and 4 columns:
#' \describe{
#'   \item{drug_name}{Generic drug name (character).}
#'   \item{brand}{Primary brand name, empty string if none noted (character).}
#'   \item{category}{Drug class (character): one of \code{"Corticosteroid"},
#'     \code{"Calcineurin inhibitor"}, \code{"PDE4 inhibitor"},
#'     \code{"JAK inhibitor"}, \code{"Biologic"}, \code{"AhR agonist"},
#'     \code{"Immunosuppressant"}.}
#'   \item{targets_known}{Semicolon-separated primary molecular target(s)
#'     from the literature (character).}
#' }
#' @source Curated from published clinical pharmacology literature (2024–2026).
#' @examples
#' data(ad_drugs)
#' ad_drugs
"ad_drugs"


#' Curated mechanism-of-action edges: AD drugs to biomarker genes
#'
#' A tibble linking each of the 18 established AD drugs to the biomarker genes
#' they regulate, with edge type and a brief mechanistic description.
#' Contains both **direct** edges (confirmed molecular target) and
#' **indirect** edges (pathway-level regulation).
#'
#' Filter to genes present in your own biomarker list before plotting:
#' ```r
#' data(ad_moa_edges)
#' ad_moa_edges |> dplyr::filter(to %in% my_gene_list)
#' ```
#'
#' @format A tibble with 106 rows and 4 columns:
#' \describe{
#'   \item{from}{Drug name — matches \code{ad_drugs$drug_name} (character).}
#'   \item{to}{HGNC gene symbol of the regulated biomarker (character).}
#'   \item{edge_type}{\code{"direct"} (confirmed molecular target) or
#'     \code{"indirect"} (pathway-level regulation).}
#'   \item{mechanism}{Brief mechanistic description of the interaction (character).}
#' }
#' @source Curated from DGIdb v5, DrugBank, and primary literature.
#' @seealso [ad_drugs] for the drug reference table.
#' @examples
#' data(ad_moa_edges)
#' ad_moa_edges |> dplyr::filter(from == "Dupilumab")
"ad_moa_edges"


#' Interaction-type colour and line-type palette
#'
#' A tibble mapping DGIdb interaction-type labels (upper-case) to hex colour
#' codes and line types used consistently across all \pkg{Repurp} plot functions.
#'
#' @format A tibble with 14 rows and 5 columns:
#' \describe{
#'   \item{interaction_type}{Interaction type label in UPPER CASE (character).}
#'   \item{edge_color}{Hex colour for network edges / links (character).}
#'   \item{sector_color}{Hex colour for chord-diagram sector band (character).}
#'   \item{link_color}{Hex colour for chord-diagram link ribbon (character).}
#'   \item{lty}{Line type integer for chord links (integer).}
#' }
#' @examples
#' data(dgi_edge_palette)
#' dgi_edge_palette
"dgi_edge_palette"


#' Drug-class colour palette for chord diagrams
#'
#' A named character vector mapping each of the 7 AD drug classes to a
#' distinct hex colour, used for sector bands in [plot_chord_dgi()].
#'
#' @format A named character vector of length 7.
#' @examples
#' data(ad_drug_class_colors)
#' ad_drug_class_colors["JAK inhibitor"]
"ad_drug_class_colors"


#' Synthetic example drug-gene interaction table
#'
#' A small, self-contained DGI data frame suitable for running examples and
#' tests without requiring external data files. Contains 15 edges linking
#' 6 AD-relevant biomarker genes to 14 drugs from both curated (MOA) and
#' DGIdb sources.
#'
#' @format A tibble with 15 rows and 4 columns:
#' \describe{
#'   \item{Gene}{HGNC gene symbol (character).}
#'   \item{Drug}{Drug name (character).}
#'   \item{InteractionType}{Interaction type string (character).}
#'   \item{Source}{\code{"Curated"} or \code{"DGIdb"} (character).}
#' }
#' @examples
#' data(example_dgi)
#' example_dgi
"example_dgi"
