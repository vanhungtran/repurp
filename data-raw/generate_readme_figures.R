# Generate README demo figures from synthetic DGI data
# Uses the same synth_dgi tribble from the vignette

suppressPackageStartupMessages({
  library(Repurp)
  library(dplyr)
  library(tibble)
  if (requireNamespace("ggplot2", quietly = TRUE)) library(ggplot2)
})

# Windows font fix
if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(sans = grDevices::windowsFont("Arial"))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::theme_set(ggplot2::theme_gray(base_family = "sans"))
  }
}

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
out <- "man/figures"

# ---- Synthetic demo data (matching vignette) ----
synth_dgi <- tribble(
  ~Gene,     ~Drug,             ~InteractionType,
  "IL13",    "Dupilumab",       "INHIBITOR",
  "IL4R",    "Dupilumab",       "INHIBITOR",
  "IL13",    "Tralokinumab",    "INHIBITOR",
  "IL13",    "Lebrikizumab",    "INHIBITOR",
  "JAK1",    "Upadacitinib",    "INHIBITOR",
  "JAK1",    "Abrocitinib",     "INHIBITOR",
  "JAK1",    "Baricitinib",     "INHIBITOR",
  "JAK2",    "Baricitinib",     "INHIBITOR",
  "JAK2",    "Ruxolitinib",     "INHIBITOR",
  "JAK1",    "Ruxolitinib",     "INHIBITOR",
  "PDE4B",   "Crisaborole",     "INHIBITOR",
  "PDE4D",   "Crisaborole",     "INHIBITOR",
  "PDE4A",   "Roflumilast",     "INHIBITOR",
  "PDE4B",   "Roflumilast",     "INHIBITOR",
  "IL31RA",  "Nemolizumab",     "INHIBITOR",
  "AHR",     "Tapinarof",       "AGONIST",
  "NR3C1",   "Hydrocortisone",  "AGONIST",
  "NR3C1",   "Triamcinolone",   "AGONIST",
  "NR3C1",   "Betamethasone",   "AGONIST",
  "FKBP1A",  "Tacrolimus",      "BINDER",
  "FKBP1A",  "Pimecrolimus",    "BINDER",
  "PPP3CA",  "Tacrolimus",      "INHIBITOR",
  "PPP3CA",  "Cyclosporine",    "INHIBITOR",
  "DHFR",    "Methotrexate",    "INHIBITOR",
  "ATIC",    "Methotrexate",    "INHIBITOR",
  "IL18",    "Hydrocortisone",  "MODULATOR",
  "IL22",    "Triamcinolone",   "MODULATOR",
  "IL33",    "Betamethasone",   "MODULATOR",
  "STAT6",   "Upadacitinib",    "INHIBITOR",
  "TSLP",    "Dupilumab",       "MODULATOR"
)

dgi_e <- transmute(synth_dgi, from = Gene, to = Drug, edge_type = InteractionType)

# ---- 1. igraph full ----
png(file.path(out, "fig_igraph_full.png"), width = 1600, height = 1400, res = 150)
plot_dgi_igraph(synth_dgi, title = "AD Drug\u2013Biomarker Network",
                gene_color = "#009E73", drug_color = "#F0E442",
                vertex_size = 7, label_cex = 0.7)
dev.off()

# ---- 2. igraph core ----
png(file.path(out, "fig_igraph_core.png"), width = 1600, height = 1400, res = 150)
plot_dgi_igraph_core(synth_dgi, title = "Core Network \u2014 genes \u2265 2 drugs",
                     min_gene_edges = 2L)
dev.off()

# ---- 3. ggraph full ----
tryCatch({
  if (requireNamespace("ggraph", quietly = TRUE)) {
    p <- plot_dgi_ggraph(dgi_e, title = "Drug\u2013Gene Network (ggraph)",
                         min_drug_targets = 1L, seed = 42L)
    ggsave(file.path(out, "fig_ggraph_full.png"), p, width = 10, height = 9, dpi = 150, bg = "white")
    message("  ggraph full: OK")
  }
}, error = function(e) message("  ggraph full: SKIP - ", e$message))

# ---- 4. ggraph enhanced (core) ----
tryCatch({
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- plot_dgi_ggraph_enhanced(dgi_e, mode = "core", title = "Enhanced Core Network",
                                  min_gene_edges = 2L, seed = 42L)
    ggsave(file.path(out, "fig_ggraph_enhanced.png"), p, width = 10, height = 9, dpi = 150, bg = "white")
    message("  ggraph enhanced: OK")
  }
}, error = function(e) message("  ggraph enhanced: SKIP - ", e$message))

# ---- 5. AD drug network ----
tryCatch({
  ad_e <- repurp_ad_edges() |>
    filter(to %in% unique(synth_dgi$Gene))

  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(ad_e) > 0) {
    p <- plot_dgi_ad_network(ad_e, drug_info = repurp_ad_drugs(),
                             title = "AD Treatment Drug\u2013Biomarker Network")
    ggsave(file.path(out, "fig_ad_network.png"), p, width = 9, height = 8, dpi = 150, bg = "white")
    message("  AD network: OK")
  }
}, error = function(e) message("  AD network: SKIP - ", e$message))

# ---- 6. Chord diagram ----
tryCatch({
  if (requireNamespace("circlize", quietly = TRUE)) {
    chord_e <- repurp_ad_edges() |>
      filter(to %in% unique(synth_dgi$Gene)) |>
      left_join(repurp_ad_drugs() |> select(drug_name, class = category),
                by = c("from" = "drug_name")) |>
      rename(drug = from, gene = to) |>
      mutate(int_type = "Curated MOA")

    svg(file.path(out, "fig_chord.svg"), width = 16, height = 16)
    plot_chord_dgi(chord_e, title = "AD Drugs \u00d7 Biomarker Chord Diagram")
    dev.off()
    message("  chord: OK")
  }
}, error = function(e) message("  chord: SKIP - ", e$message))

# ---- 7. Heatmap ----
tryCatch({
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- plot_dgi_heatmap(dgi_e, drug_info = repurp_ad_drugs(),
                          title = "Drug \u00d7 Gene Interaction Heatmap",
                          max_drugs = 20L, max_genes = 15L)
    ggsave(file.path(out, "fig_heatmap.png"), p, width = 11, height = 8, dpi = 150, bg = "white")
    message("  heatmap: OK")
  }
}, error = function(e) message("  heatmap: SKIP - ", e$message))

# ---- 8. Dot matrix ----
tryCatch({
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- plot_dgi_dotmatrix(dgi_e, drug_info = repurp_ad_drugs(),
                            title = "Gene \u00d7 Drug-Class Connectivity", max_genes = 15L)
    ggsave(file.path(out, "fig_dotmatrix.png"), p, width = 9, height = 8, dpi = 150, bg = "white")
    message("  dot matrix: OK")
  }
}, error = function(e) message("  dot matrix: SKIP - ", e$message))

# ---- 9. Sankey ----
tryCatch({
  if (requireNamespace("ggalluvial", quietly = TRUE)) {
    p <- plot_dgi_sankey(dgi_e, drug_info = repurp_ad_drugs(),
                         fill_by = "drug_class", title = "Drug Class \u2192 Drug \u2192 Gene Flow",
                         max_genes = 15L)
    ggsave(file.path(out, "fig_sankey.png"), p, width = 10, height = 8, dpi = 150, bg = "white")
    message("  sankey: OK")
  }
}, error = function(e) message("  sankey: SKIP - ", e$message))

# ---- 10. Circular arc ----
tryCatch({
  if (requireNamespace("ggraph", quietly = TRUE)) {
    circ_e <- repurp_ad_edges() |>
      filter(to %in% unique(synth_dgi$Gene)) |>
      left_join(repurp_ad_drugs() |> select(drug_name, class = category),
                by = c("from" = "drug_name")) |>
      rename(class = category)

    p <- plot_dgi_circular_arc(circ_e, title = "AD Drug\u2013Biomarker Circular Network")
    ggsave(file.path(out, "fig_circular_arc.png"), p, width = 9, height = 9, dpi = 150, bg = "white")
    message("  circular arc: OK")
  }
}, error = function(e) message("  circular arc: SKIP - ", e$message))

cat("\nFigures saved to man/figures/\n")
list.files(out)
