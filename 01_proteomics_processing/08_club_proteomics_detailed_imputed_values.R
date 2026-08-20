library(writexl)

ctrl <- read.csv("secretory_cntrl_b17_b20.csv")
il13 <- read.csv("secretory_il13_f2_f5.csv")

analyze_protein <- function(ctrl_vals, il13_vals) {
  
  # ----------------------------
  # Step 1: Log2 transform
  # ----------------------------
  control_log <- log2(ctrl_vals)
  treated_log <- log2(il13_vals)
  
  # keep track of which values were originally missing
  ctrl_was_imputed <- is.na(control_log)
  il13_was_imputed <- is.na(treated_log)
  
  # ----------------------------
  # Step 2: Compute imputation value
  # ----------------------------
  observed_values <- c(control_log, treated_log)
  observed_values <- observed_values[!is.na(observed_values)]
  
  if (length(observed_values) < 2) {
    return(c(
      Log2.Mean.ctrl = NA,
      Log2.Mean.il13 = NA,
      Log2.Difference = NA,
      p.value = NA,
      ctrl_1 = NA, ctrl_2 = NA, ctrl_3 = NA, ctrl_4 = NA,
      il13_1 = NA, il13_2 = NA, il13_3 = NA, il13_4 = NA,
      ctrl_imp_1 = NA, ctrl_imp_2 = NA, ctrl_imp_3 = NA, ctrl_imp_4 = NA,
      il13_imp_1 = NA, il13_imp_2 = NA, il13_imp_3 = NA, il13_imp_4 = NA
    ))
  }
  
  mean_obs <- mean(observed_values)
  sd_obs <- sd(observed_values)
  impute_value <- mean_obs - 1.8 * sd_obs
  
  # ----------------------------
  # Step 3: Impute missing values
  # ----------------------------
  if (any(is.na(control_log))) {
    control_log[is.na(control_log)] <- rnorm(sum(is.na(control_log)), mean = impute_value, sd = 0.1)
  }
  
  if (any(is.na(treated_log))) {
    treated_log[is.na(treated_log)] <- rnorm(sum(is.na(treated_log)), mean = impute_value, sd = 0.1)
  }
  
  # ----------------------------
  # Step 4: Summary stats
  # ----------------------------
  mean_ctrl <- mean(control_log)
  mean_il13 <- mean(treated_log)
  p_value <- tryCatch(t.test(control_log, treated_log)$p.value, error = function(e) NA)
  
  # pad to length 4 in case some sheets have only 3 replicates
  pad_to_4 <- function(x) {
    length(x) <- 4
    x
  }
  
  control_log <- pad_to_4(control_log)
  treated_log <- pad_to_4(treated_log)
  ctrl_was_imputed <- pad_to_4(ctrl_was_imputed)
  il13_was_imputed <- pad_to_4(il13_was_imputed)
  
  # return results
  return(c(
    Log2.Mean.ctrl = mean_ctrl,
    Log2.Mean.il13 = mean_il13,
    Log2.Difference = mean_il13 - mean_ctrl,
    p.value = p_value,
    ctrl_1 = control_log[1],
    ctrl_2 = control_log[2],
    ctrl_3 = control_log[3],
    ctrl_4 = control_log[4],
    il13_1 = treated_log[1],
    il13_2 = treated_log[2],
    il13_3 = treated_log[3],
    il13_4 = treated_log[4],
    ctrl_imp_1 = ctrl_was_imputed[1],
    ctrl_imp_2 = ctrl_was_imputed[2],
    ctrl_imp_3 = ctrl_was_imputed[3],
    ctrl_imp_4 = ctrl_was_imputed[4],
    il13_imp_1 = il13_was_imputed[1],
    il13_imp_2 = il13_was_imputed[2],
    il13_imp_3 = il13_was_imputed[3],
    il13_imp_4 = il13_was_imputed[4]
  ))
}

set.seed(123)

results_mtrx <- matrix(NA, nrow = nrow(ctrl), ncol = 20)

colnames(results_mtrx) <- c(
  "Log2.Mean.ctrl",
  "Log2.Mean.il13",
  "Log2.Difference",
  "p.value",
  "ctrl_1",
  "ctrl_2",
  "ctrl_3",
  "ctrl_4",
  "il13_1",
  "il13_2",
  "il13_3",
  "il13_4",
  "ctrl_imp_1",
  "ctrl_imp_2",
  "ctrl_imp_3",
  "ctrl_imp_4",
  "il13_imp_1",
  "il13_imp_2",
  "il13_imp_3",
  "il13_imp_4"
)

for (x in 1:nrow(ctrl)) {
  ctrl_vals <- as.numeric(ctrl[x, 5:8])
  il13_vals <- as.numeric(il13[x, 5:8])
  
  results_mtrx[x, ] <- analyze_protein(ctrl_vals, il13_vals)
}

results_detailed_df <- as.data.frame(results_mtrx)

# add back protein data

results_detailed_df$Protein.Group <- ctrl$Protein.Group
results_detailed_df$Protein.Names <- ctrl$Protein.Names
results_detailed_df$Genes <- ctrl$Genes
results_detailed_df$Description <- ctrl$First.Protein.Description


results_detailed_df <- results_detailed_df[, c(
  "Protein.Group",
  "Protein.Names",
  "Genes",
  "Description",
  "Log2.Mean.ctrl",
  "Log2.Mean.il13",
  "Log2.Difference",
  "p.value",
  "ctrl_1",
  "ctrl_2",
  "ctrl_3",
  "ctrl_4",
  "il13_1",
  "il13_2",
  "il13_3",
  "il13_4",
  "ctrl_imp_1",
  "ctrl_imp_2",
  "ctrl_imp_3",
  "ctrl_imp_4",
  "il13_imp_1",
  "il13_imp_2",
  "il13_imp_3",
  "il13_imp_4"
)]

write_xlsx(results_detailed_df, "secretory_imputed_detailed_results.xlsx")
write.csv(results_detailed_df, "secretory_imputed_detailed_results.csv", row.names = FALSE)
