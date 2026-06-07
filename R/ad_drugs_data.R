#' Curated reference table of 18 established atopic-dermatitis (AD) drugs
#'
#' Returns a tibble with drug name, brand name, drug class, known molecular
#' targets, and mechanism-of-action notes.
#'
#' @return A tibble with columns `drug_name`, `brand`, `category`, `targets_known`.
#' @export
#' @examples
#' repurp_ad_drugs()
repurp_ad_drugs <- function() {
  tibble::tribble(
    ~drug_name,       ~brand,       ~category,                    ~targets_known,
    "Hydrocortisone", "",           "Corticosteroid",             "NR3C1 (GR)",
    "Triamcinolone",  "",           "Corticosteroid",             "NR3C1 (GR)",
    "Betamethasone",  "",           "Corticosteroid",             "NR3C1 (GR)",
    "Tacrolimus",     "Protopic",   "Calcineurin inhibitor",      "FKBP1A; PPP3CA",
    "Pimecrolimus",   "Elidel",     "Calcineurin inhibitor",      "FKBP1A; PPP3CA",
    "Cyclosporine",   "",           "Calcineurin inhibitor",      "PPP3CA; PPP3CB",
    "Crisaborole",    "Eucrisa",    "PDE4 inhibitor",             "PDE4B; PDE4D",
    "Roflumilast",    "Zoryve",     "PDE4 inhibitor",             "PDE4A; PDE4B; PDE4D",
    "Ruxolitinib",    "Opzelura",   "JAK inhibitor",              "JAK1; JAK2",
    "Upadacitinib",   "Rinvoq",     "JAK inhibitor",              "JAK1",
    "Abrocitinib",    "Cibinqo",    "JAK inhibitor",              "JAK1",
    "Baricitinib",    "Olumiant",   "JAK inhibitor",              "JAK1; JAK2",
    "Dupilumab",      "Dupixent",   "Biologic",                   "IL4R",
    "Tralokinumab",   "Adbry",      "Biologic",                   "IL13",
    "Lebrikizumab",   "Ebglyss",    "Biologic",                   "IL13",
    "Nemolizumab",    "Nemluvio",   "Biologic",                   "IL31RA",
    "Tapinarof",      "Vtama",      "AhR agonist",                "AHR",
    "Methotrexate",   "",           "Immunosuppressant",          "DHFR; ATIC"
  )
}


#' Curated mechanism-of-action edges: AD drugs → biomarker genes
#'
#' Returns a long-format tibble linking each of the 18 established AD drugs to
#' the biomarker genes they regulate, with edge type (`"direct"` or `"indirect"`)
#' and a brief mechanism string.
#'
#' Filter to the genes present in your own biomarker list before plotting:
#' ```r
#' repurp_ad_edges() |> dplyr::filter(to %in% my_gene_list)
#' ```
#'
#' @return A tibble with columns `from` (drug name), `to` (gene symbol),
#'   `edge_type` (`"direct"` / `"indirect"`), and `mechanism`.
#' @export
#' @examples
#' repurp_ad_edges()
repurp_ad_edges <- function() {
  dplyr::bind_rows(
    # ---- Corticosteroids ----
    tibble::tibble(from = "Hydrocortisone",
                   to   = c("IL13","IL22","IL18","IL33","IL4R","IL2RA",
                             "IFN-gamma","CCL22","CCL28","CCL23","MCP-4"),
                   edge_type = "indirect",
                   mechanism = "GR -> down-NF-kB -> down-cytokine"),
    tibble::tibble(from = "Triamcinolone",
                   to   = c("IL13","IL22","IL18","IL33","IL4R",
                             "IL2RA","IFN-gamma","MCP-4"),
                   edge_type = "indirect",
                   mechanism = "GR -> down-NF-kB"),
    tibble::tibble(from = "Betamethasone",
                   to   = c("IL13","IL22","IL18","IL33","IL4R",
                             "IL2RA","IFN-gamma","MCP-4"),
                   edge_type = "indirect",
                   mechanism = "GR -> down-NF-kB"),

    # ---- Calcineurin inhibitors ----
    tibble::tibble(from = "Tacrolimus",
                   to   = c("IL13","IL22","IL4R","IL2RA","IFN-gamma","IL18","CCL22"),
                   edge_type = "indirect",
                   mechanism = "down-Calcineurin -> down-NFAT -> down-cytokine"),
    tibble::tibble(from = "Pimecrolimus",
                   to   = c("IL13","IL22","IL4R","IL2RA","IFN-gamma"),
                   edge_type = "indirect",
                   mechanism = "down-Calcineurin -> down-NFAT"),
    tibble::tibble(from = "Cyclosporine",
                   to   = c("IL13","IL22","IL4R","IL2RA","IFN-gamma","IL18","CCL22"),
                   edge_type = "indirect",
                   mechanism = "down-Calcineurin -> down-NFAT -> down-cytokine"),

    # ---- PDE4 inhibitors ----
    tibble::tibble(from = "Crisaborole",
                   to   = c("IL13","IL22","IL4R","IFN-gamma","IL18"),
                   edge_type = "indirect",
                   mechanism = "up-cAMP -> down-cytokine"),
    tibble::tibble(from = "Roflumilast",
                   to   = c("IL13","IL22","IL4R","IFN-gamma","IL18","IL33","CCL22"),
                   edge_type = "indirect",
                   mechanism = "up-cAMP -> down-cytokine"),

    # ---- JAK inhibitors ----
    tibble::tibble(from = "Ruxolitinib",
                   to   = c("IL4R","IL13","IL22","IL2RA","IFN-gamma","IL18","IL33"),
                   edge_type = "indirect",
                   mechanism = "down-JAK1/2 -> down-STAT -> down-cytokine signaling"),
    tibble::tibble(from = "Upadacitinib",
                   to   = c("IL4R","IL13","IL22","IL2RA","IFN-gamma","IL18","IL33"),
                   edge_type = "indirect",
                   mechanism = "down-JAK1 -> down-STAT -> down-cytokine signaling"),
    tibble::tibble(from = "Abrocitinib",
                   to   = c("IL4R","IL13","IL22","IL2RA","IFN-gamma","IL18","IL33"),
                   edge_type = "indirect",
                   mechanism = "down-JAK1 -> down-STAT -> down-cytokine signaling"),
    tibble::tibble(from = "Baricitinib",
                   to   = c("IL4R","IL13","IL22","IL2RA","IFN-gamma","IL18","IL33"),
                   edge_type = "indirect",
                   mechanism = "down-JAK1/2 -> down-STAT -> down-cytokine signaling"),

    # ---- Biologics (direct molecular targets) ----
    tibble::tibble(from = "Dupilumab",    to = "IL4R",  edge_type = "direct",
                   mechanism = "anti-IL4R-alpha mAb"),
    tibble::tibble(from = "Tralokinumab", to = "IL13",  edge_type = "direct",
                   mechanism = "anti-IL13 mAb"),
    tibble::tibble(from = "Lebrikizumab", to = "IL13",  edge_type = "direct",
                   mechanism = "anti-IL13 mAb"),
    tibble::tibble(from = "Nemolizumab",
                   to   = c("IL4R","IL13","IL22"),
                   edge_type = "indirect",
                   mechanism = "anti-IL31RA -> down-type-2 cross-talk"),

    # ---- AhR agonist ----
    tibble::tibble(from = "Tapinarof",
                   to   = c("IL22","IL13","IL4R","IL18"),
                   edge_type = "indirect",
                   mechanism = "AhR -> down-IL-22/Th17"),

    # ---- Immunosuppressant ----
    tibble::tibble(from = "Methotrexate",
                   to   = c("IL13","IL22","IL4R","IL2RA","IFN-gamma",
                             "IL18","IL33","CCL22","CCL28","MCP-4"),
                   edge_type = "indirect",
                   mechanism = "up-adenosine -> down-inflammation")
  )
}
