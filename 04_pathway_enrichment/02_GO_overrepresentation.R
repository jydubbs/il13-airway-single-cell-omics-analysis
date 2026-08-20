# For figures 4 and 5 - Category-specific pathway enrichment (GO:BP)

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# figure4_genesets is the folder that has all of my gene sets for each category prepared
setwd("figure4_genesets")

run_go <- function(foreground_file, background_file) {
  
  foreground <- unique(na.omit(read.csv(foreground_file)$gene))
  background <- unique(na.omit(read.csv(background_file)$gene))
  foreground <- intersect(foreground, background)
  
  ego <- enrichGO(
    gene          = foreground,
    universe      = background,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10,
    minGSSize     = 5,
    maxGSSize     = 500
  )
  
  # Stop if no GO terms were enriched
  if (nrow(as.data.frame(ego)) == 0) {
    return(NULL)
  }
  
  # Remove highly similar/redundant GO terms
  ego <- simplify(
    ego,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  return(ego)
}

# ------------------------------------------------------------
# The three cell types. Each has its own background file:
#   bg_goblet.csv  bg_club.csv  bg_ciliated.csv
# ------------------------------------------------------------
cells      <- c("goblet", "club", "ciliated")

# for plotting
nice_names <- c("Goblet", "Club", "Multiciliated")

# ------------------------------------------------------------
# Loops over the three cell types, runs run_go() for each,
# then puts them side by side with merge_result() + dotplot().
# ------------------------------------------------------------
make_panel <- function(list_prefix, out_png, out_csv, plot_title) {
  results <- list()
  # Run ORA separately for each cell type
  for (i in seq_along(cells)) {
    cell <- cells[i]
    nicename <- nice_names[i]
    fg_file <- paste0(list_prefix, "_", cell, ".csv")
    bg_file <- paste0("bg4_", cell, ".csv")
    ego <- run_go(fg_file, bg_file)
    # Only keep cell types where enrichment was found
    if (!is.null(ego)) {
      results[[nicename]] <- ego
    }
  }
  
  # Combine cell types into one comparison plot
  cc <- merge_result(results)
  # Save enrichment results
  write.csv(
    as.data.frame(cc),
    out_csv,
    row.names = FALSE
  )
  
  # Plot
  p <- dotplot(cc, showCategory = 8) +
    ggtitle(plot_title)
  print(p)
  ggsave(
    out_png,
    p,
    width = 9,
    height = 10,
    dpi = 300
  )

}

# ============================================================
# Build each panel by calling make_panel()
# ============================================================

# --- 4A: Concordant, UP with IL13 ---
make_panel("fg4_concordant_up",
           "fig4A_concordant_up.png",
           "fig4A_concordant_up_GOBP.csv",
           "Concordant UP with IL13")

# --- 4B: Concordant, DOWN with IL13 ---
make_panel("fg4_concordant_down",
           "fig4B_concordant_down.png",
           "fig4B_concordant_down_GOBP.csv",
           "Concordant DOWN with IL13")

# --- 5A: Protein-only, UP with IL13 ---
make_panel("fg5_protein_only_up",
           "fig5A_protein_only_up.png",
           "fig5A_protein_only_up_GOBP.csv",
           "Protein-only UP with IL13")

# --- 5B: Protein-only, DOWN with IL13 ---
make_panel("fg5_protein_only_down",
           "fig5B_protein_only_down.png",
           "fig5B_protein_only_down_GOBP.csv",
           "Protein-only DOWN with IL13")

