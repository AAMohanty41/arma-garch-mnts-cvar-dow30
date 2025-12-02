#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
})

message("=== 10_backtest_plots_2025.R : plotting 2025 backtest results ===")

# --------------------------------------------------------------------
# 0) Load data
# --------------------------------------------------------------------
sum_path <- "outputs/tables/backtest_2025_summary.csv"
win_path <- "outputs/tables/backtest_2025_windows_detail.csv"

if (!file.exists(sum_path) || !file.exists(win_path)) {
  stop("Missing backtest CSVs. Run 09_compare_2025.R first.")
}

sum25 <- read_csv(sum_path, show_col_types = FALSE)
win25 <- read_csv(win_path, show_col_types = FALSE)

sum25$portfolio <- factor(
  sum25$portfolio,
  levels = c("MS_historical", "MS_CVaR_ARMA", "MS_CVaR_noARMA")
)

dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------------------------
# 1) Predicted vs realized mean 21d return
# --------------------------------------------------------------------
mean_long <- sum25 %>%
  select(portfolio, pred_mean_21d, realized_mean_21d) %>%
  pivot_longer(
    cols = c(pred_mean_21d, realized_mean_21d),
    names_to = "type",
    values_to = "mean_21d"
  ) %>%
  mutate(type = recode(type,
                       pred_mean_21d     = "Predicted",
                       realized_mean_21d = "Realized"))

p_mean <- ggplot(mean_long, aes(x = portfolio, y = mean_21d, fill = type)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Predicted vs realized 21-day mean return (2025)",
    x = NULL,
    y = "21-day log return"
  ) +
  theme_minimal()

ggsave("outputs/figures/backtest2025_pred_vs_real_mean.png",
       p_mean, width = 7, height = 4, dpi = 300)

# --------------------------------------------------------------------
# 2) Predicted vs realized 21d volatility (SD)
# --------------------------------------------------------------------
sd_long <- sum25 %>%
  select(portfolio, pred_sd_21d, realized_sd_21d) %>%
  pivot_longer(
    cols = c(pred_sd_21d, realized_sd_21d),
    names_to = "type",
    values_to = "sd_21d"
  ) %>%
  mutate(type = recode(type,
                       pred_sd_21d     = "Predicted",
                       realized_sd_21d = "Realized"))

p_sd <- ggplot(sd_long, aes(x = portfolio, y = sd_21d, fill = type)) +
  geom_col(position = "dodge") +
  labs(
    title = "Predicted vs realized 21-day volatility (2025)",
    x = NULL,
    y = "SD of 21-day log return"
  ) +
  theme_minimal()

ggsave("outputs/figures/backtest2025_pred_vs_real_sd.png",
       p_sd, width = 7, height = 4, dpi = 300)

# --------------------------------------------------------------------
# 3) VaR breach rate (how calibrated are the models?)
# --------------------------------------------------------------------
p_var <- ggplot(sum25, aes(x = portfolio, y = VaR95_breach_rate)) +
  geom_col() +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(
    title = "Empirical 95% VaR breach rate in 2025",
    x = NULL,
    y = "Fraction of 21-day windows with loss > VaR"
  ) +
  theme_minimal()

ggsave("outputs/figures/backtest2025_var_breach_rate.png",
       p_var, width = 7, height = 4, dpi = 300)

# --------------------------------------------------------------------
# 4) Time series of realized 21d returns across 2025
# --------------------------------------------------------------------
win_long <- win25 %>%
  mutate(
    start_date = as.Date(start_date),
    end_date   = as.Date(end_date),
    mid_date   = start_date + (end_date - start_date) / 2
  ) %>%
  select(window_id, mid_date,
         realized_ms_21d,
         realized_cvar_arma_21d,
         realized_cvar_noarma_21d) %>%
  pivot_longer(
    cols = starts_with("realized_"),
    names_to = "portfolio",
    values_to = "ret_21d"
  ) %>%
  mutate(portfolio = recode(
    portfolio,
    realized_ms_21d          = "MS_historical",
    realized_cvar_arma_21d   = "MS_CVaR_ARMA",
    realized_cvar_noarma_21d = "MS_CVaR_noARMA"
  ))

p_ts <- ggplot(win_long, aes(x = mid_date, y = ret_21d, color = portfolio)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Realized 21-day portfolio returns across 2025",
    x = "Window midpoint",
    y = "21-day log return"
  ) +
  theme_minimal()

ggsave("outputs/figures/backtest2025_realized_returns_timeseries.png",
       p_ts, width = 7.5, height = 4.5, dpi = 300)

# --------------------------------------------------------------------
# 5) Percentile distribution (calibration of full MNTS sims)
# --------------------------------------------------------------------
perc_long <- win25 %>%
  select(window_id,
         percentile_ms,
         percentile_cvar_arma,
         percentile_cvar_noarma) %>%
  pivot_longer(
    cols = -window_id,
    names_to = "portfolio",
    values_to = "percentile"
  ) %>%
  mutate(portfolio = recode(
    portfolio,
    percentile_ms          = "MS_historical",
    percentile_cvar_arma   = "MS_CVaR_ARMA",
    percentile_cvar_noarma = "MS_CVaR_noARMA"
  ))

p_perc <- ggplot(perc_long, aes(x = portfolio, y = percentile)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(
    title = "Distribution of realized 21-day percentiles vs MNTS sims (2025)",
    x = NULL,
    y = "Percentile in simulated distribution"
  ) +
  theme_minimal()

ggsave("outputs/figures/backtest2025_percentile_boxplot.png",
       p_perc, width = 7, height = 4, dpi = 300)

message("=== Done: 10_backtest_plots_2025.R (figures in outputs/figures/) ===")
