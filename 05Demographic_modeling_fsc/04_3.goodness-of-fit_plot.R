#!/usr/bin/env Rscript
# =============================================================================
# Goodness of Fit Test for fastsimcoal2 Model Validation
# =============================================================================

library(ggplot2)
library(reshape2)
library(gridExtra)

# =============================================================================
# 1. Read Data
# =============================================================================
setwd("E://paper/07 costatae evolution/data_analysis/fastsimcoal 20250409/Goodness_of_fit")
# Read simulated data (5000 rows)
sim_data <- read.table("all.stats.modelFixT-031.sim.txt", header = TRUE, sep = "\t")
obs_data <- read.table("observed_data.txt", header = TRUE, sep = "\t")


# =============================================================================
# 2. Calculate P-values (Quantile Position)
# =============================================================================

# Function: Calculate two-tailed p-value
calc_pvalue <- function(obs, sim_vector) {
  n <- length(sim_vector)
  rank_lower <- sum(sim_vector <= obs) / n
  rank_upper <- sum(sim_vector >= obs) / n
  
  # Two-tailed p-value
  p_value <- 2 * min(rank_lower, rank_upper)
  p_value <- min(p_value, 1)  # Ensure it does not exceed 1
  
  return(list(
    quantile = rank_lower,
    p_value = p_value
  ))
}

# Calculate p-value for each statistic
results <- data.frame(
  Statistic = colnames(sim_data),
  Observed = numeric(ncol(sim_data)),
  Sim_Mean = numeric(ncol(sim_data)),
  Sim_SD = numeric(ncol(sim_data)),
  Sim_Min = numeric(ncol(sim_data)),
  Sim_Max = numeric(ncol(sim_data)),
  Quantile = numeric(ncol(sim_data)),
  P_value = numeric(ncol(sim_data)),
  Significant = character(ncol(sim_data))
)

for (i in 1:ncol(sim_data)) {
  stat_name <- colnames(sim_data)[i]
  obs_val <- obs_data[[stat_name]]
  sim_vec <- sim_data[[stat_name]]
  
  pval_result <- calc_pvalue(obs_val, sim_vec)
  
  results$Observed[i] <- obs_val
  results$Sim_Mean[i] <- mean(sim_vec)
  results$Sim_SD[i] <- sd(sim_vec)
  results$Sim_Min[i] <- min(sim_vec)
  results$Sim_Max[i] <- max(sim_vec)
  results$Quantile[i] <- pval_result$quantile
  results$P_value[i] <- pval_result$p_value
  
  # Mark significance
  if (pval_result$p_value < 0.01) {
    results$Significant[i] <- "***"
  } else if (pval_result$p_value < 0.05) {
    results$Significant[i] <- "**"
  } else if (pval_result$p_value < 0.1) {
    results$Significant[i] <- "*"
  } else {
    results$Significant[i] <- ""
  }
}


# =============================================================================
# 3. Plot Histograms
# =============================================================================

plot_histogram <- function(stat_name, sim_vec, obs_val, pval) {
  df <- data.frame(value = sim_vec)
  
  p <- ggplot(df, aes(x = value)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
    geom_vline(xintercept = obs_val, color = "red", linewidth = 1.2, linetype = "dashed") +
    labs(
      title = stat_name,
      subtitle = sprintf("p = %.4f", pval),
      x = "Value",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8, color = ifelse(pval < 0.05, "red", "black")),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7)
    )
  
  return(p)
}

# Plot histograms for all statistics
plot_list <- list()
for (i in 1:ncol(sim_data)) {
  stat_name <- colnames(sim_data)[i]
  plot_list[[i]] <- plot_histogram(
    stat_name,
    sim_data[[stat_name]],
    obs_data[[stat_name]],
    results$P_value[i]
  )
}

# Plot Pi and FST separately
# Pi statistics
pi_idx <- grep("^Pi_", colnames(sim_data))
if (length(pi_idx) > 0) {
  pi_plots <- plot_list[pi_idx]
  n_col <- min(3, length(pi_plots))
  n_row <- ceiling(length(pi_plots) / n_col)
  
  pdf("goodness_of_fit_Pi.pdf", width = 10, height = 3 * n_row)
  do.call(grid.arrange, c(pi_plots, ncol = n_col))
  dev.off()
  cat("\nPi statistics plot saved to: goodness_of_fit_Pi.pdf\n")
}

# FST statistics
fst_idx <- grep("^FST_", colnames(sim_data))
if (length(fst_idx) > 0) {
  fst_plots <- plot_list[fst_idx]
  n_col <- 3
  n_row <- ceiling(length(fst_plots) / n_col)
  
  pdf("goodness_of_fit_FST.pdf", width = 12, height = 3 * n_row)
  do.call(grid.arrange, c(fst_plots, ncol = n_col))
  dev.off()
  cat("FST statistics plot saved to: goodness_of_fit_FST.pdf\n")
}

# Combined plot for all statistics
pdf("goodness_of_fit_all.pdf", width = 15, height = 15)
n_col <- 4
n_row <- ceiling(length(plot_list) / n_col)
do.call(grid.arrange, c(plot_list, ncol = n_col))
dev.off()
cat("All statistics plot saved to: goodness_of_fit_all.pdf\n")

