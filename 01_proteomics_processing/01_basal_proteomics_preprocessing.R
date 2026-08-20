library(writexl)

ctrl <- read.csv("basal_cntrl_b7_b10.csv")
il13 <- read.csv("basal_il13_f12_f15.csv")

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
# Step 9: Adding missingness classes
# ---------------------------- 

sig_hits_annotated$class <- paste(sig_hits_annotated$observed_ctrl, sig_hits_annotated$observed_il13, sep="-")
sig_hits_annotated <- sig_hits_annotated[order(sig_hits_annotated$p.value), ]

# ---------------------------- 
# Step 10: Export as xls and csv files
# ---------------------------- 

write_xlsx(sig_hits_annotated, "ver_3_basal_sig_hits_results.xlsx")
write.csv(sig_hits_annotated, "ver_3_basal_sig_hits_results.csv", row.names = FALSE)

