
ctrl <- read.csv(path_1)
il13 <- read.csv(path_2)

# first define function for imputing values

analyze_protein <- function(ctrl_vals, il13_vals) {
  
  # ---------------------------- 
  # Step 1: Log2 transform 
  # ---------------------------- 
  control_log <- log2(ctrl_vals)
  treated_log <- log2(il13_vals)
  
  # ---------------------------- 
  # Step 2: Compute imputation value 
  # ----------------------------
  
  # Take all non-NA observed values from both groups 
  observed_values <- c(control_log, treated_log)
  observed_values <- observed_values[!is.na(observed_values)]
  
  # if too few values to get mean and sd (e.g. 1 or 0 total control+il13 values), return NA for everything
  if (length(observed_values) < 2) {
    return(c(mean_ctrl = NA, mean_il13 = NA, log2FC = NA, pvalue = NA))
  }
  
  mean_obs <- mean(observed_values)
  sd_obs <- sd(observed_values)
  
  # imputation value
  impute_value <- mean_obs - 1.8 * sd_obs
  
  # ---------------------------- 
  # Step 3: Impute missing values 
  # ---------------------------- 
  
  # to generate exactly as many random numbers (around the imputation value) as there are missing values in results (using normal distribution)
  if (any(is.na(control_log))) {
    control_log[is.na(control_log)] <- rnorm(sum(is.na(control_log)), mean = impute_value, sd = 0.1)
  }
  
  if (any(is.na(treated_log))) {
    treated_log[is.na(treated_log)] <- rnorm(sum(is.na(treated_log)), mean = impute_value, sd = 0.1)
  }
  
  # get log2 means
  mean_ctrl <- mean(control_log)
  mean_il13 <- mean(treated_log)
  
  # t-test
  p_value <- tryCatch(t.test(control_log, treated_log)$p.value, error = function(e) NA)
  
  # return results
  return(c(Log2.Mean.ctrl = mean_ctrl,
           Log2.Mean.il13 = mean_il13,
           Log2.Difference = mean_il13 - mean_ctrl,
           p.value = p_value))
}

# ---------------------------- 
# Step 4: Pipeline for imputation + analysis (log2 means, log2 fold change, t-test)
# ---------------------------- 

set.seed(123)

results_mtrx <- matrix(NA, nrow = nrow(ctrl), ncol = 4)

colnames(results_mtrx) <- c(
  "Log2.Mean.ctrl",
  "Log2.Mean.il13",
  "Log2.Difference",
  "p.value"
)

for (x in 1:nrow(ctrl)) {
  ctrl_vals <- as.numeric(ctrl[x, 5:8])
  il13_vals <- as.numeric(il13[x, 5:8])
  
  results_mtrx[x, ] <- analyze_protein(ctrl_vals, il13_vals)
}

results_df <- as.data.frame(results_mtrx)

# ---------------------------- 
# Step 5: Add in other metadata
# ---------------------------- 

# corresponding protein info
results_df$Protein.Group <- ctrl$Protein.Group
results_df$Protein.Names <- ctrl$Protein.Names
results_df$Genes <- ctrl$Genes
results_df$Description <- ctrl$First.Protein.Description

# obtain ranking by significance (p.value)
results_df <- results_df[order(results_df$p.value), ]
results_df$Rank <- seq_len(nrow(results_df))

# assuming false discovery rate of 0.05 - get the BH critical value
results_df$BH.critical.value <- (results_df$Rank / nrow(results_df)) * 0.05

# if the p value is less than BH critical value, call significance = "Significant", if not, NA
results_df$Significance <- ifelse(
  results_df$p.value <= results_df$BH.critical.value, "Significant", "Not significant")

results_df$Significant.up..1 <- ifelse(
  results_df$Significance == "Significant" & results_df$Log2.Difference >= 1,
  results_df$Genes, NA)
results_df$Significant.down..1 <- ifelse(
  results_df$Significance == "Significant" & results_df$Log2.Difference <= -1,
  results_df$Genes, NA)
results_df$Significant.up..0.5 <- ifelse(
  results_df$Significance == "Significant" & results_df$Log2.Difference >= 0.5,
  results_df$Genes, NA)
results_df$Significant.down..0.5 <- ifelse(
  results_df$Significance == "Significant" & results_df$Log2.Difference <= -0.5,
  results_df$Genes, NA)

# mean abundance to find treand in protein expression levels
results_df$Log2.Mean <- (results_df$Log2.Mean.ctrl + results_df$Log2.Mean.il13) / 2

# reorder everything to look cleaner
final_df <- results_df[, c(
  "Protein.Group",
  "Protein.Names",
  "Genes",
  "Description",
  "Log2.Mean.ctrl",
  "Log2.Mean.il13",
  "Log2.Mean",
  "Log2.Difference",
  "p.value",
  "Rank",
  "BH.critical.value",
  "Significance",
  "Significant.up..1",
  "Significant.down..1",
  "Significant.up..0.5",
  "Significant.down..0.5"
)]

# ---------------------------- 
# Step 6: Import data sets
# ---------------------------- 

sig_hits <- final_df

# check that the rows for the ctrl and il13 are in the same order
# extract first column of each matrix
ctrl_1 <- ctrl[, 1]
il13_1 <- il13[, 1]

# compare the sorted columns
result <- identical(ctrl_1, il13_1)
print(result) # should be TRUE

# ---------------------------- 
# Step 7: Deciding which proteins/groups have low observed values (1 or less values for each group)
# ---------------------------- 

# select the columns with data
ctrl_values <- ctrl[,5:8]
il13_values <- il13[,5:8]

# converts values to TRUE/FALSE depending on observed counts, then count number of observed values
observed_ctrl <- rowSums(!is.na(ctrl_values))
observed_il13 <- rowSums(!is.na(il13_values))

# identify any columns with only 1 or less observed values
flag_ctrl_high_missing <- observed_ctrl <= 1
flag_il13_high_missing <- observed_il13 <= 1

# if there the are 0 values for ctrl and 2 or more values for il13 = flag control sided absence (may be imputation-driven results, but biologically true)- then vice verca
flag_ctrl_sided_absent <- observed_ctrl == 0 & observed_il13 >= 2
flag_il13_sided_absent <- observed_il13 == 0 & observed_ctrl >= 2

# ---------------------------- 
# Step 8: Concatenating to sig_hits dataframe (with all the analysis) by Protein.Group
# ---------------------------- 

# make new dataframe - on Protein.Group
flags_df <- data.frame(
  Protein.Group = ctrl$Protein.Group,
  observed_ctrl = observed_ctrl,
  observed_il13 = observed_il13,
  flag_ctrl_high_missing = flag_ctrl_high_missing,
  flag_il13_high_missing = flag_il13_high_missing,
  flag_ctrl_sided_absent = flag_ctrl_sided_absent,
  flag_il13_sided_absent = flag_il13_sided_absent,
  stringsAsFactors = FALSE
)

# merge tables - attach metadata to significant proteins by ID
sig_hits_annotated <- merge(
  sig_hits,
  flags_df,
  by = "Protein.Group",
  all.x = TRUE
)

# ---------------------------- 
# Step 9: Adding classes
# ---------------------------- 

sig_hits_annotated$class <- paste(sig_hits_annotated$observed_ctrl, sig_hits_annotated$observed_il13, sep="-")

#interesting_classes <- c("2-2", "0-4", "4-0", "3-3", "1-4", "4-1", "4-4")
interesting_classes <- c("2-2", "0-4", "4-0", "3-3", "1-4", "4-1", "4-4")

plot_df <- sig_hits_annotated[
  sig_hits_annotated$class %in% interesting_classes,
]

plot_df$class <- factor(
  plot_df$class,
  levels = c("2-2", "0-4", "4-0", "3-3", "1-4", "4-1", "4-4")
)

library(ggplot2)

ggplot(plot_df, aes(x = class, y = Log2.Mean)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue") +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(
    title = "Protein abundance distribution by missingness class
                          in Basal Stem Cells",
    x = "Class (ctrl-il13)",
    y = "Log2 Mean Abundance"
  ) +
  theme_minimal()

# Splitting up and down regulated

plot_df$direction <- ifelse(plot_df$Log2.Difference > 0, "Up", "Down")

ggplot(plot_df, aes(x = class, y = Log2.Mean, fill = direction)) +
  geom_boxplot() +
  theme_minimal()

# Conclusion:
# 4-4 proteins occupy the highest abundance tier:
# - dense distribution + median around 17-18 + wide dynamic range up to 30
# - lots of overlap between upregulated and downregulated proteins and both high abundance
# - IL-13 is regulating already abundant proteins (core biology)
# other classes:
# - lower/narrower spread + lower medians + lower upper range + range between 12-23 + tight around ~12-15
# - 4-1 proteins are significantly lower abundance → supports detection-driven missingness
# - not much trend between upregulated and downregulated proteins
# - the separation graph shows generally lower abundance, slight but weak separatation between up vs down
# - 4-1/1-4 likely borderline detection proteins, more esensitive to imputation strategies -> leading to significance
# 1-4 proteins tend to occupy lower abundance ranges but partially overlap with core biological pathways
# identified in 4-4, suggesting they may represent lower-abundance components of the same biological
# programs rather than entirely distinct processes.

# adding stats - shows the visual
kruskal.test(Log2.Mean ~ class, data = plot_df)
pairwise.wilcox.test(plot_df$Log2.Mean, plot_df$class)

# show variance in the different classes - should see lower variance in imputation heavy-classes
tapply(plot_df$Log2.Mean, plot_df$class, sd)

theme_map <- list(
  lipid_metabolism = c("FA2H","ELOVL5","FUT6","NEU3","SCARB1","STARD4","ATP10B","PCTP","ACSL3","CPT1A"),
  membrane_trafficking = c("SFTPD","DLL1","VAMP5","VAMP8","SNAP23","STX7","RAB11FIP5","EHD1","NSF","SCAMP4"),
  cytoskeleton_adhesion = c("CSPG4","AIF1L","CDC42EP5","FLRT3","SEMA3F","ITGAV","ITGB5","VCL","FLNB","ROCK2"),
  immune_signaling = c("SFTPD","HLA-DRB5","IL18","TLR3","TAP1","IFIT1","IRF9","NFKBIA","ICAM1","CD44"),
  apoptosis_mito = c("BNIP3","BAX","CASP3","CASP7","BID","AIFM1","CYCS","TXN","GPX4","SOD1"),
  protein_modification = c("ZDHHC8","VCPKMT","P3H2","PRMT1","PRMT5","HDAC6","DNMT1","UBE2G1","OTUB1","USP28"),
  nucleotide_metabolism = c("NUDT3","NUDT4B","NUDT4","NT5C","PNP","IMPDH2","DHFR","GART","PAICS"),
  cell_cycle = c("E2F6","CDK1","CDK4","PCNA","MCM2","MCM5","MCM7","BUB3","RFC4"),
  differentiation = c("DLL1","HOXA1","RUNX2","SOX2","KRT5","KRT13","IVL","DSP","PKP1"),
  ion_signaling = c("ITPR1","CAMKK1","STIM1","TRPM4","KCNQ1","PIEZO1"),
  metabolism_general = c("CSAD","GGACT","GSTA2","GCLC","MAT2A","TPI1","LDHA","HK1"),
  unknown = c("BAALC","MOB3B","ZSCAN18","ZNF649","ZBTB48","MANSC1","YAE1")
)

sig_hits_annotated$theme <- "unknown"

assign_theme <- function(gene_string, theme_map) {
  gene_vec <- unlist(strsplit(as.character(gene_string), ";"))
  gene_vec <- trimws(gene_vec)
  
  matched_themes <- names(theme_map)[sapply(theme_map, function(x) any(gene_vec %in% x))]
  
  if (length(matched_themes) == 0) {
    return("unknown")
  } else {
    return(matched_themes[1])
  }
}

sig_hits_annotated$theme <- sapply(
  sig_hits_annotated$Genes,
  assign_theme,
  theme_map = theme_map
)

plot_df <- sig_hits_annotated[
  sig_hits_annotated$class %in% interesting_classes,
]

plot_df$class <- factor(
  plot_df$class,
  levels = c("2-2", "0-4", "4-0", "3-3", "1-4", "4-1", "4-4"),
  ordered = TRUE
)


# table 1: class vs abundance

sig_hits_annotated$abundance_bin <- cut(
  sig_hits_annotated$Log2.Mean,
  breaks = quantile(sig_hits_annotated$Log2.Mean,
                    probs = c(0, 0.25, 0.5, 0.75, 1),
                    na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("Low", "Mid-low", "Mid-high", "High")
)

plot_df <- sig_hits_annotated[
  sig_hits_annotated$class %in% interesting_classes,
]

#  Now run table
class_abundance <- prop.table(
  table(plot_df$class, plot_df$abundance_bin),
  1
)
write.csv(as.data.frame.matrix(class_abundance),
          "/Users/janicewong/Documents/University/Queen Mary/Research project/Data/sc Proteomics/Analysis/basal/basal_class_abundance.csv")


# table 2: class vs theme
class_theme_prop <- prop.table(table(plot_df$class, plot_df$theme), 1)
class_theme_counts <- table(plot_df$class, plot_df$theme)

write.csv(as.data.frame.matrix(class_theme_prop),
          "/Users/janicewong/Documents/University/Queen Mary/Research project/Data/sc Proteomics/Analysis/basal/basal_class_theme_prop.csv")

write.csv(as.data.frame.matrix(class_theme_counts),
          "/Users/janicewong/Documents/University/Queen Mary/Research project/Data/sc Proteomics/Analysis/basal/basal_class_theme_counts.csv")


# table 3: theme vs abundance

theme_abundance_prop <- prop.table(table(plot_df$theme, plot_df$abundance_bin), 1)
theme_abundance_counts <- table(plot_df$theme, plot_df$abundance_bin)

write.csv(as.data.frame.matrix(theme_abundance_prop),
          "/Users/janicewong/Documents/University/Queen Mary/Research project/Data/sc Proteomics/Analysis/basal/basal_theme_abundance_prop.csv")

write.csv(as.data.frame.matrix(theme_abundance_counts),
          "/Users/janicewong/Documents/University/Queen Mary/Research project/Data/sc Proteomics/Analysis/basal/basal_theme_abundance_counts.csv")

#### additional plots

ggplot(sig_hits_annotated, aes(x = theme, y = Log2.Mean)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_theme_df <- sig_hits_annotated[
  sig_hits_annotated$class %in% c("2-2","3-3","0-4","4-0","1-4","4-1","4-4") &
    sig_hits_annotated$theme != "unknown",
]

plot_theme_df$class <- factor(
  plot_theme_df$class,
  levels = c("2-2","3-3","0-4","4-0","1-4","4-1","4-4")
)

ggplot(plot_theme_df, aes(x = class, fill = theme)) +
  geom_bar(position = "fill") +
  theme_minimal()

ggplot(plot_theme_df, aes(x = class, fill = theme)) +
  geom_bar(position = "stack") +
  theme_minimal()

# ---------------------------- 
# Step 9b: Add 3-level abundance bins
# ---------------------------- 

# binning according to protein expression levels
make_bins <- function(x) {
  cut(x,
    breaks = quantile(x,probs = c(0, 1/3, 2/3, 1),na.rm = TRUE),
    include.lowest = TRUE,
    labels = c("Low", "Mid", "High"))
}

sig_hits_annotated$ctrl_abundance_bin <- make_bins(sig_hits_annotated$Log2.Mean.ctrl)
sig_hits_annotated$il13_abundance_bin <- make_bins(sig_hits_annotated$Log2.Mean.il13)
sig_hits_annotated$mean_abundance_bin <- make_bins(sig_hits_annotated$Log2.Mean)

# ---------------------------- 
# Step 10: Interpreting the flags
# ---------------------------- 

#sig_hits_annotated$interpretation <- "High confidence"

# moderate confidence
# if either condition is poorly observed = interpret cautiously e.g. control: NA NA value, il13: value value value
#sig_hits_annotated$interpretation[
#  (sig_hits_annotated$flag_ctrl_high_missing | sig_hits_annotated$flag_il13_high_missing)
#]  <- "Moderate confidence"

# low confidence
# if both conditions are porrly observed = very low confidence e.g. control: NA NA value, il13: value NA NA
#sig_hits_annotated$interpretation[
#  sig_hits_annotated$flag_ctrl_high_missing & sig_hits_annotated$flag_il13_high_missing
#]  <- "Low confidence (high missingness)"

# is there a one-sided absence? could be imputation-driven significance (p-value), but biologically true
#sig_hits_annotated$interpretation[
#  sig_hits_annotated$flag_ctrl_sided_absent
#] <- "Significance may be imputation-driven (absent control values)"

#sig_hits_annotated$interpretation[
#  sig_hits_annotated$flag_il13_sided_absent
#] <- "Significance may be imputation-driven (absent IL13-treated values)"

sig_hits_annotated <- sig_hits_annotated[order(sig_hits_annotated$p.value), ]



library(writexl)

write_xlsx(sig_hits_annotated, "path_3")
write.csv(sig_hits_annotated, "path_4", row.names = FALSE)

write_xlsx(sig_hits_annotated, "path_5")
write.csv(sig_hits_annotated, "path_6", row.names = FALSE)

