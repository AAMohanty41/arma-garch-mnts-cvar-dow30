#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

message("=== 11_summarize_2025_sharpe_drawdown.R ===")

summary_path <- "outputs/tables/backtest_2025_summary.csv"
windows_path <- "outputs/tables/backtest_2025_windows_detail.csv"

if (!file.exists(summary_path)) {
  stop("Missing ", summary_path, " – run the 2025 backtest scripts first.")
}
if (!file.exists(windows_path)) {
  stop("Missing ", windows_path, " – run the 2025 backtest scripts first.")
}

summary_2025 <- read_csv(summary_path, show_col_types = FALSE)
windows_2025 <- read_csv(windows_path, show_col_types = FALSE)

# -----------------------------
# 1) Realized Sharpe (annualized)
# -----------------------------
h_days        <- 21          # horizon
trading_days  <- 252
scale_mu      <- trading_days / h_days
scale_sd      <- sqrt(trading_days / h_days)

# assume rf ~ 0 for 2025
summary_2025 <- summary_2025 %>%
  mutate(
    realized_sharpe_annual =
      (realized_mean_21d * scale_mu) /
      (realized_sd_21d   * scale_sd)
  )

sharpe_out <- summary_2025 %>%
  select(
    portfolio,
    realized_mean_21d,
    realized_sd_21d,
    realized_sharpe_annual
  )

write_csv(
  sharpe_out,
  "outputs/tables/backtest_2025_realized_sharpe.csv"
)

message("Wrote realized Sharpe table -> outputs/tables/backtest_2025_realized_sharpe.csv")

# -----------------------------
# 2) Max drawdown from 21d windows
# -----------------------------

compute_max_dd <- function(ret_21d_vec) {
  wealth <- exp(cumsum(ret_21d_vec))   # turn 21d log returns into wealth path
  peak   <- cummax(wealth)
  dd     <- wealth / peak - 1          # drawdown series
  min(dd)                              # most negative drawdown
}

dd_tbl <- tibble(
  portfolio    = c("MS_historical", "MS_CVaR_ARMA", "MS_CVaR_noARMA"),
  max_drawdown = c(
    compute_max_dd(windows_2025$realized_ms_21d),
    compute_max_dd(windows_2025$realized_cvar_arma_21d),
    compute_max_dd(windows_2025$realized_cvar_noarma_21d)
  )
)

write_csv(
  dd_tbl,
  "outputs/tables/backtest_2025_max_drawdown.csv"
)

message("Wrote max drawdown table -> outputs/tables/backtest_2025_max_drawdown.csv")
message("=== done 11_summarize_2025_sharpe_drawdown.R ===")
