.repurp_global_variables <- c(
  ".data", "Drug", "Gene", "InteractionType", "brand", "category",
  "color", "coverage", "deg", "degree", "directional_edges", "drug",
  "drug_class", "drug_count", "drug_name", "edge_alpha", "edge_category",
  "edge_count", "edge_linetype", "edge_reversal_score", "edge_score",
  "edge_type", "edge_width", "effect_direction", "fdr", "freq", "from",
  "from_id", "from_name", "gene", "gene_score", "group", "id", "idx",
  "int_type", "interaction_type", "interactions", "is_hub", "label",
  "label2", "label_face", "label_size", "n", "n_connections", "n_drugs",
  "name", "node_fill", "point_size", "reversal_score", "score", "sources",
  "start", "stratum", "supportive_edges", "target_count", "targets", "to",
  "to_id", "to_name", "total", "type", "width", "x", "x_from", "x_strip",
  "x_to", "y", "y_from", "y_to"
)

utils::globalVariables(.repurp_global_variables)

.onLoad <- function(libname, pkgname) {
  if (.Platform$OS.type == "windows") {
    grDevices::windowsFonts(sans = grDevices::windowsFont("Arial"))
  }
}
