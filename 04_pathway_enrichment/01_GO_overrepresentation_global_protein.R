# ============================================================
# Figure 3 - GO:BP enrichment of protein-significant genes
# Split by IL13 direction:
#   Figure 3A = genes UP with IL13
#   Figure 3B = genes DOWN with IL13
# Goblet, Club, Multiciliated compared side by side in each panel.
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

setwd("figure_3_genesets")

# ------------------------------------------------------------
# Small helper: run GO:BP enrichment for ONE cell type against
# ITS OWN background, then merge near-duplicate terms.
# foreground_file = protein-significant genes for that cell type + direction
# background_file = all protein-measured genes for that cell type
# ------------------------------------------------------------
run_go <- function(foreground_file, background_file) {
  foreground <- read.csv(foreground_file)$gene
  background <- read.csv(background_file)$gene
  ego <- enrichGO(
    gene         = foreground,
    OrgDb        = org.Hs.eg.db,
    keyType      = "SYMBOL",
    ont          = "BP",
    universe     = background,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.10
  )
  simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
}

# ============================================================
# Figure 3A - genes UP with IL13
# ============================================================
up_goblet <- run_go("fg_up_goblet.csv",   "bg_goblet.csv")
up_club   <- run_go("fg_up_club.csv",     "bg_club.csv")
up_multi  <- run_go("fg_up_ciliated.csv", "bg_ciliated.csv")

# Combine the three cell types into one side-by-side result
cc_up <- merge_result(list(
  Goblet        = up_goblet,
  Club          = up_club,
  Multiciliated = up_multi
))

write.csv(as.data.frame(cc_up), "fig3a_up_GOBP.csv", row.names = FALSE)

p_up <- dotplot(cc_up, showCategory = 8) +
  ggtitle("Figure 3A - GO:BP enriched in proteins UP with IL13")
ggsave("fig3a_up.png", p_up, width = 9, height = 10, dpi = 300)

print(p_up)

# ============================================================
# Figure 3B - genes DOWN with IL13
# ============================================================
down_goblet <- run_go("fg_down_goblet.csv",   "bg_goblet.csv")
down_club   <- run_go("fg_down_club.csv",     "bg_club.csv")
down_multi  <- run_go("fg_down_ciliated.csv", "bg_ciliated.csv")

cc_down <- merge_result(list(
  Goblet        = down_goblet,
  Club          = down_club,
  Multiciliated = down_multi
))

write.csv(as.data.frame(cc_down), "fig3b_down_GOBP.csv", row.names = FALSE)

p_down <- dotplot(cc_down, showCategory = 8) +
  ggtitle("Figure 3B - GO:BP enriched in proteins DOWN with IL13")
ggsave("fig3b_down.png", p_down, width = 9, height = 10, dpi = 300)

print(p_down)
