#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

message("=== 08_plot_cvar_frontier.R : plotting CVaR frontier ===")

# Read the CVaR frontier table we just created
frontier_cvar <- read.csv("outputs/tables/frontier_cvar_21d.csv",
                          stringsAsFactors = FALSE)

if (!all(c("mean_21d", "sd_21d", "CVaR95_loss", "Sharpe_sd") %in% names(frontier_cvar))) {
  stop("frontier_cvar_21d.csv is missing required columns. Check 07_frontier_cvar.R output.")
}

# Identify max-Sharpe (sd-based) point on the CVaR frontier
best_idx <- which.max(frontier_cvar$Sharpe_sd)

frontier_cvar$point_type <- "frontier"
if (length(best_idx) == 1 && is.finite(frontier_cvar$Sharpe_sd[best_idx])) {
  frontier_cvar$point_type[best_idx] <- "max_sharpe"
} else {
  warning("Could not cleanly identify a unique max-Sharpe point; plotting all points as 'frontier'.")
}

# -------------------------------------------------------------------
# Plot 1: mean vs CVaR (this is the “mean–CVaR frontier”)
# -------------------------------------------------------------------

p_cvar <- ggplot(frontier_cvar,
                 aes(x = CVaR95_loss, y = mean_21d)) +
  geom_line(alpha = 0.6) +
  geom_point(aes(shape = point_type), size = 2.3) +
  scale_shape_manual(values = c(frontier = 16, max_sharpe = 17)) +
  labs(
    title = "21-day Mean–CVaR Frontier (MNTS, ARMA–GARCH)",
    x = "CVaR 95% (loss, 21-day)",
    y = "Expected return (21-day)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

ggsave("outputs/figures/frontier_cvar_21d.png",
       plot   = p_cvar,
       width  = 7,
       height = 5,
       dpi    = 300)

message("Saved CVaR frontier figure to outputs/figures/frontier_cvar_21d.png")

message("=== Done: 08_plot_cvar_frontier.R ===")
