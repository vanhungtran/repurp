# Run once to (re-)generate all .rda files under data/
# Rscript data-raw/build_package_data.R

library(tibble)
library(dplyr)

# ============================================================
#  1. ad_drugs  —  18 curated AD drugs reference table
# ============================================================
ad_drugs <- tribble(
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

# ============================================================
#  2. ad_moa_edges  —  curated MOA edges (drug -> gene)
# ============================================================
ad_moa_edges <- bind_rows(
  # Corticosteroids
  tibble(from = "Hydrocortisone",
         to   = c("IL13","IL22","IL18","IL33","IL4R","IL2RA",
                  "IFN-gamma","CCL22","CCL28","CCL23","MCP-4"),
         edge_type = "indirect", mechanism = "GR -> down-NF-kB -> down-cytokine"),
  tibble(from = "Triamcinolone",
         to   = c("IL13","IL22","IL18","IL33","IL4R","IL2RA","IFN-gamma","MCP-4"),
         edge_type = "indirect", mechanism = "GR -> down-NF-kB"),
  tibble(from = "Betamethasone",
         to   = c("IL13","IL22","IL18","IL33","IL4R","IL2RA","IFN-gamma","MCP-4"),
         edge_type = "indirect", mechanism = "GR -> down-NF-kB"),
  # Calcineurin inhibitors
  tibble(from = "Tacrolimus",
         to   = c("IL13","IL22","IL4R","IL2RA","IFN-gamma","IL18","CCL22"),
         edge_type = "indirect", mechanism = "down-Calcineurin -> down-NFAT -> down-cytokine"),
  tibble(from = "Pimecrolimus",
         to   = c("IL13","IL22","IL4R","IL2RA","IFN-gamma"),
         edge_type = "indirect", mechanism = "down-Calcineurin -> down-NFAT"),
  tibble(from = "Cyclosporine",
         to   = c("IL13","IL22","IL4R","IL2RA","IFN-gamma","IL18","CCL22"),
         edge_type = "indirect", mechanism = "down-Calcineurin -> down-NFAT -> down-cytokine"),
  # PDE4 inhibitors
  tibble(from = "Crisaborole",
         to   = c("IL13","IL22","IL4R","IFN-gamma","IL18"),
         edge_type = "indirect", mechanism = "up-cAMP -> down-cytokine"),
  tibble(from = "Roflumilast",
         to   = c("IL13","IL22","IL4R","IFN-gamma","IL18","IL33","CCL22"),
         edge_type = "indirect", mechanism = "up-cAMP -> down-cytokine"),
  # JAK inhibitors
  tibble(from = "Ruxolitinib",
         to   = c("IL4R","IL13","IL22","IL2RA","IFN-gamma","IL18","IL33"),
         edge_type = "indirect", mechanism = "down-JAK1/2 -> down-STAT -> down-cytokine"),
  tibble(from = "Upadacitinib",
         to   = c("IL4R","IL13","IL22","IL2RA","IFN-gamma","IL18","IL33"),
         edge_type = "indirect", mechanism = "down-JAK1 -> down-STAT -> down-cytokine"),
  tibble(from = "Abrocitinib",
         to   = c("IL4R","IL13","IL22","IL2RA","IFN-gamma","IL18","IL33"),
         edge_type = "indirect", mechanism = "down-JAK1 -> down-STAT -> down-cytokine"),
  tibble(from = "Baricitinib",
         to   = c("IL4R","IL13","IL22","IL2RA","IFN-gamma","IL18","IL33"),
         edge_type = "indirect", mechanism = "down-JAK1/2 -> down-STAT -> down-cytokine"),
  # Biologics
  tibble(from = "Dupilumab",    to = "IL4R",
         edge_type = "direct",   mechanism = "anti-IL4R-alpha mAb"),
  tibble(from = "Tralokinumab", to = "IL13",
         edge_type = "direct",   mechanism = "anti-IL13 mAb"),
  tibble(from = "Lebrikizumab", to = "IL13",
         edge_type = "direct",   mechanism = "anti-IL13 mAb"),
  tibble(from = "Nemolizumab",
         to   = c("IL4R","IL13","IL22"),
         edge_type = "indirect", mechanism = "anti-IL31RA -> down-type-2 cross-talk"),
  # AhR agonist
  tibble(from = "Tapinarof",
         to   = c("IL22","IL13","IL4R","IL18"),
         edge_type = "indirect", mechanism = "AhR -> down-IL-22/Th17"),
  # Immunosuppressant
  tibble(from = "Methotrexate",
         to   = c("IL13","IL22","IL4R","IL2RA","IFN-gamma",
                  "IL18","IL33","CCL22","CCL28","MCP-4"),
         edge_type = "indirect", mechanism = "up-adenosine -> down-inflammation")
)

# ============================================================
#  3. dgi_edge_palette  — interaction-type colour/linetype table
# ============================================================
dgi_edge_palette <- tribble(
  ~interaction_type,  ~edge_color,  ~sector_color, ~link_color, ~lty,
  "INHIBITOR",        "#E41A1C",    "#FB9A99",    "#FBB4AE",   2L,
  "AGONIST",          "#377EB8",    "#A6CEE3",    "#B3CDE3",   3L,
  "ANTAGONIST",       "#FF7F00",    "#FDBF6F",    "#FED9A6",   4L,
  "ACTIVATOR",        "#4DAF4A",    "#B2DF8A",    "#E5F5E0",   2L,
  "BINDER",           "#984EA3",    "#CAB2D6",    "#DECBE4",   5L,
  "BLOCKER",          "#A65628",    "#B15928",    "#CCEBC5",   6L,
  "MODULATOR",        "#F781BF",    "#FCCDE5",    "#FFF7BC",   3L,
  "SUBSTRATE",        "#66C2A5",    "#80CDC1",    "#C2E699",   4L,
  "ANTIBODY",         "#8DD3C7",    "#8DD3C7",    "#A6D854",   5L,
  "CLEAVAGE",         "#BEBADA",    "#BEBADA",    "#D4B9DA",   6L,
  "INVERSE AGONIST",  "#FB8072",    "#FB8072",    "#FDD0A2",   2L,
  "OTHER",            "#BDBDBD",    "#BDBDBD",    "#D9D9D9",   1L,
  "NONE",             "#999999",    "#CCCCCC",    "#DDDDDD",   1L,
  "CURATED MOA",      "#404040",    "#404040",    "#888888",   1L
)

# ============================================================
#  4. ad_drug_class_colors  — drug-class sector colour map
# ============================================================
ad_drug_class_colors <- c(
  "Corticosteroid"         = "#E41A1C",
  "Calcineurin inhibitor"  = "#FF7F00",
  "PDE4 inhibitor"         = "#377EB8",
  "JAK inhibitor"          = "#4DAF4A",
  "Biologic"               = "#984EA3",
  "AhR agonist"            = "#F781BF",
  "Immunosuppressant"      = "#A65628"
)

# ============================================================
#  5. example_dgi  —  small synthetic DGI table for examples/tests
# ============================================================
example_dgi <- tibble(
  Gene            = c("IL13","IL13","IL13","IL4R","IL4R","IL22","IL22",
                      "IL18","IL18","IL33","NTRK1","NTRK1","EGF","EGF","CCL22"),
  Drug            = c("Dupilumab","Tralokinumab","Lebrikizumab",
                      "Dupilumab","Baricitinib","Tapinarof","Roflumilast",
                      "Baricitinib","Ruxolitinib","Upadacitinib",
                      "Imatinib","Entrectinib","Cetuximab","Erlotinib",
                      "Methotrexate"),
  InteractionType = c("direct","direct","direct",
                      "direct","INHIBITOR","indirect","indirect",
                      "INHIBITOR","INHIBITOR","INHIBITOR",
                      "INHIBITOR","INHIBITOR","ANTAGONIST","INHIBITOR",
                      "indirect"),
  Source          = c(rep("Curated", 10), rep("DGIdb", 5))
)

# ============================================================
#  Save all objects
# ============================================================
save(ad_drugs,           file = "data/ad_drugs.rda",           compress = "xz")
save(ad_moa_edges,       file = "data/ad_moa_edges.rda",       compress = "xz")
save(dgi_edge_palette,   file = "data/dgi_edge_palette.rda",   compress = "xz")
save(ad_drug_class_colors, file = "data/ad_drug_class_colors.rda", compress = "xz")
save(example_dgi,        file = "data/example_dgi.rda",        compress = "xz")

cat("Package data written to data/\n")
cat("  ad_drugs            ", nrow(ad_drugs), "rows\n")
cat("  ad_moa_edges        ", nrow(ad_moa_edges), "rows\n")
cat("  dgi_edge_palette    ", nrow(dgi_edge_palette), "rows\n")
cat("  ad_drug_class_colors", length(ad_drug_class_colors), "entries\n")
cat("  example_dgi         ", nrow(example_dgi), "rows\n")
