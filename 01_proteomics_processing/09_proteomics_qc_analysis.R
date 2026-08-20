library(ggplot2)

# ============================================================
# 1. Load data
# ============================================================

goblet_ctrl <- read.csv("goblet_cntrl_b12_b15.csv")
goblet_il13 <- read.csv("goblet_il13_f7_f10.csv")

basal_ctrl <- read.csv("basal_cntrl_b7_b10.csv")
basal_il13 <- read.csv("basal_il13_f12_f15.csv")

mccs_ctrl <- read.csv("mccs_ctrl_b2_b5.csv")
mccs_il13 <- read.csv("mccs_il-13_f17_f20.csv")

secretory_ctrl <- read.csv("secretory_cntrl_b17_b20.csv")
secretory_il13 <- read.csv("secretory_il13_f2_f5.csv")

muco_cil_ctrl <- read.csv("muco_cil_cntr_d2_d5.csv")
muco_cil_il13 <- read.csv("muco_cil_il13_h2_h5.csv")

muco_basal_ctrl <- read.csv("muco_basal_cntrl_d7_d10.csv")
muco_basal_il13 <- read.csv("muco_basal_il13_h7_h10.csv")


# ============================================================
# 2. Flexible replicate extraction
# ============================================================

extract_replicates_flexible <- function(df, celltype, condition) {
  expr <- df[, 5:ncol(df)]
  rownames(expr) <- df$Protein.Group
  n_reps <- ncol(expr)
  colnames(expr) <- paste(
    celltype,
    condition,
    paste0("rep", 1:n_reps),
    sep = "_"
  )
  return(expr)
}


# ============================================================
# 3. Extract matrices
# ============================================================

goblet_ctrl_expr <- extract_replicates_flexible(goblet_ctrl, "goblet", "control")
goblet_il13_expr <- extract_replicates_flexible(goblet_il13, "goblet", "IL13")

basal_ctrl_expr <- extract_replicates_flexible(basal_ctrl, "basal", "control")
basal_il13_expr <- extract_replicates_flexible(basal_il13, "basal", "IL13")

mccs_ctrl_expr <- extract_replicates_flexible(mccs_ctrl, "mccs", "control")
mccs_il13_expr <- extract_replicates_flexible(mccs_il13, "mccs", "IL13")

secretory_ctrl_expr <- extract_replicates_flexible(secretory_ctrl, "secretory", "control")
secretory_il13_expr <- extract_replicates_flexible(secretory_il13, "secretory", "IL13")

muco_cil_ctrl_expr <- extract_replicates_flexible(muco_cil_ctrl, "muco_cil", "control")
muco_cil_il13_expr <- extract_replicates_flexible(muco_cil_il13, "muco_cil", "IL13")

muco_basal_ctrl_expr <- extract_replicates_flexible(muco_basal_ctrl, "muco_basal", "control")
muco_basal_il13_expr <- extract_replicates_flexible(muco_basal_il13, "muco_basal", "IL13")


# ============================================================
# 4. PCA helper function
# ============================================================

run_proteomics_pca <- function(expr_list, title_text, keep_threshold = 0.7) {

  all_expr <- do.call(cbind, expr_list)
  all_expr <- as.matrix(all_expr)
  mode(all_expr) <- "numeric"

  all_log <- log2(all_expr)
  all_log[!is.finite(all_log)] <- NA

  keep <- rowMeans(!is.na(all_log)) >= keep_threshold
  all_log_filt <- all_log[keep, ]

  cat("\n", title_text, "\n")
  cat("Proteins before filtering:", nrow(all_log), "\n")
  cat("Proteins used for PCA:", nrow(all_log_filt), "\n")

  set.seed(123)
  all_log_imp <- all_log_filt

  for (i in 1:nrow(all_log_imp)) {

    obs <- all_log_imp[i, !is.na(all_log_imp[i, ])]

    if (length(obs) >= 2) {

      mean_obs <- mean(obs)
      sd_obs <- sd(obs)

      impute_value <- mean_obs - 1.8 * sd_obs
      n_missing <- sum(is.na(all_log_imp[i, ]))

      if (n_missing > 0) {
        all_log_imp[i, is.na(all_log_imp[i, ])] <- rnorm(
          n_missing,
          mean = impute_value,
          sd = 0.1
        )
      }
    }
  }

  all_log_imp <- all_log_imp[
    complete.cases(all_log_imp) &
      apply(all_log_imp, 1, function(x) all(is.finite(x))),
  ]

  pca <- prcomp(
    t(all_log_imp),
    center = TRUE,
    scale. = TRUE
  )

  pca_df <- data.frame(
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )

  sample_split <- strsplit(pca_df$Sample, "_")

  pca_df$Celltype <- sapply(sample_split, function(x) {
    paste(x[1:(length(x) - 2)], collapse = "_")
  })

  pca_df$Condition <- sapply(sample_split, function(x) {
    x[length(x) - 1]
  })

  pca_df$Replicate <- sapply(sample_split, function(x) {
    x[length(x)]
  })

  pc_var <- summary(pca)$importance[2, 1:2] * 100

  p <- ggplot(
    pca_df,
    aes(
      x = PC1,
      y = PC2,
      colour = Celltype,
      shape = Condition
    )
  ) +
    geom_point(size = 4, alpha = 0.9) +
    theme_minimal(base_size = 14) +
    labs(
      title = title_text,
      x = paste0("PC1 (", round(pc_var[1], 1), "%)"),
      y = paste0("PC2 (", round(pc_var[2], 1), "%)")
    )

  return(list(
    pca = pca,
    pca_df = pca_df,
    plot = p
  ))
}


# ============================================================
# 5. PCA for 4 major epithelial populations only
# ============================================================

major_pca <- run_proteomics_pca(
  expr_list = list(
    goblet_ctrl_expr,
    goblet_il13_expr,
    basal_ctrl_expr,
    basal_il13_expr,
    mccs_ctrl_expr,
    mccs_il13_expr,
    secretory_ctrl_expr,
    secretory_il13_expr
  ),
  title_text = "PCA of proteomics samples: major epithelial populations",
  keep_threshold = 0.7
)

major_pca$plot +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "right"
  ) +
  coord_cartesian(xlim = c(-90, 90))

# ============================================================
# 6. Protein groups quantified per well (depth QC)
# ============================================================
# NOTE: use the RAW _expr matrices here (non-NA counts), NOT the
# imputed/filtered PCA matrix — depth = what was actually measured.

count_depth <- function(expr, celltype, condition) {
  data.frame(
    Sample     = colnames(expr),
    Celltype   = celltype,
    Condition  = condition,
    n_proteins = colSums(!is.na(expr))
  )
}

depth_all <- rbind(
  count_depth(goblet_ctrl_expr,    "goblet",    "control"),
  count_depth(goblet_il13_expr,    "goblet",    "IL13"),
  count_depth(secretory_ctrl_expr, "secretory", "control"),
  count_depth(secretory_il13_expr, "secretory", "IL13"),
  count_depth(mccs_ctrl_expr,      "mccs",      "control"),
  count_depth(mccs_il13_expr,      "mccs",      "IL13"),
  count_depth(basal_ctrl_expr,     "basal",     "control"),
  count_depth(basal_il13_expr,     "basal",     "IL13")
)

# fix the order cell types appear on the x-axis, and condition colour order
depth_all$Celltype  <- factor(depth_all$Celltype,
                              levels = c("goblet", "secretory", "mccs", "basal"))
depth_all$Condition <- factor(depth_all$Condition,
                              levels = c("control", "IL13"))

# ============================================================
# 7. Plot: bars = group mean, points = individual wells
# ============================================================

# grand mean across all wells — compute AFTER depth_all has all cell types rbind-ed
grand_mean <- mean(depth_all$n_proteins)

p_depth <- ggplot(depth_all, aes(x = Celltype, y = n_proteins, fill = Condition)) +
  stat_summary(fun = mean, geom = "col",
               position = position_dodge(width = 0.8),
               width = 0.7, alpha = 0.55) +
  geom_point(aes(colour = Condition),
             position = position_dodge(width = 0.8),
             size = 1.0,
             shape  = 21,
             colour = "black",
             stroke = 0.4,
             show.legend = FALSE) +
  geom_hline(yintercept = grand_mean, linetype = "dashed",
             colour = "grey30", linewidth = 0.6) +
  annotate("text", x = 3.8, y = grand_mean,
           label = paste0("mean = ", round(grand_mean)),
           hjust = 0, vjust = -0.5, size = 3.5, colour = "grey30") +
  scale_fill_manual(values  = c(control = "#3E8E5A", IL13 = "#E1812C")) +
  scale_colour_manual(values = c(control = "#3E8E5A", IL13 = "#E1812C")) +
  scale_y_continuous(breaks       = seq(0, 10000, 2000),
                     minor_breaks = seq(0, 10000, 1000),
                     limits       = c(0, 10000),
                     expand       = expansion(mult = c(0, 0.02))) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Protein groups per well",
    x = NULL,
    y = "No. protein groups per well"
  ) +
  theme(
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    legend.position    = "right",
    panel.grid.minor.y = element_line(linewidth = 0.3, colour = "grey88"),
    panel.grid.major.y = element_line(linewidth = 0.5, colour = "grey80")
  )

p_depth

# ============================================================
# 8. Purity - marker genes heatmap
# ============================================================

xl_path <- "report.pg_matrix_pplate3_groups.xlsx"

goblet   <- read_excel(xl_path, "GOBLET B12_B15")
club     <- read_excel(xl_path, "SECRETORY B17_B20")
ciliated <- read_excel(xl_path, "MCCS CTR B2_B5")
basal    <- read_excel(xl_path, "BASAL B7_B10")

markers <- c("MUC5AC","MUC5B","FCGBP","BPIFB1", "TSPAN8",     # goblet
             "SCGB1A1","BPIFA1","CEACAM5","CEACAM6",  # club
             "FOXJ1","DNAI1","DNAI2","RSPH1",         # ciliated
             "KRT5","KRT14","TP63","NGFR")            # basal

get_means <- function(df) {
  intensity <- log2(as.matrix(df[, 5:ncol(df)]))
  intensity[!is.finite(intensity)] <- NA
  means <- rowMeans(intensity, na.rm = TRUE)
  means[match(markers, df$Genes)]
}

mat <- data.frame(Goblet = get_means(goblet), Club = get_means(club),
                  Ciliated = get_means(ciliated), Basal = get_means(basal))
rownames(mat) <- markers

z <- as.data.frame(t(scale(t(mat))))   # z-score each row
z$Marker <- markers

long <- pivot_longer(z, cols = c("Goblet","Club","Ciliated","Basal"),
                     names_to = "Population", values_to = "zscore")
long$Marker     <- factor(long$Marker, levels = rev(markers))
long$Population  <- factor(long$Population, levels = c("Goblet","Club","Ciliated","Basal"))

ggplot(long, aes(Population, Marker, fill = zscore)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#4A6FA5", mid = "white", high = "#B0413E", midpoint = 0, name = "z-score") +
  labs(x = "", y = "", fill = "z-score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

