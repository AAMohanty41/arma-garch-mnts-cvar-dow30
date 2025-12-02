#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(tidyquant)
  library(dplyr)
  library(tidyr)
  library(readr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

message("=== 09_compare_2025.R : compare predicted vs realized 21d windows in 2025 ===")

# ------------------------------------------------------------
# 1) Load config + MNTS sims
# ------------------------------------------------------------
cfg     <- yaml::read_yaml("config/config.yml")
tickers <- cfg$universe
alpha   <- cfg$alpha %||% 0.95

sim_path_arma   <- "data/processed/sims_mnts_21d.rds"
sim_path_noarma <- "data/processed/sims_mnts_21d_noarma.rds"

if (!file.exists(sim_path_arma) || !file.exists(sim_path_noarma)) {
  stop("Missing sims_mnts_21d*.rds — run 06_mnts_fit_sim.R and 06b_mnts_fit_sim_noarma.R first.")
}

R_sims_arma   <- as.matrix(readRDS(sim_path_arma))
R_sims_noarma <- as.matrix(readRDS(sim_path_noarma))

S        <- nrow(R_sims_arma)
n_assets <- ncol(R_sims_arma)

sim_tickers <- colnames(R_sims_arma)
if (is.null(sim_tickers)) {
  stop("sims_mnts_21d.rds has no column names.")
}

if (!setequal(sim_tickers, tickers)) {
  warning("Tickers in sims vs cfg$universe differ; using sims' tickers.")
}
tickers       <- sim_tickers
R_sims_arma   <- R_sims_arma[,   tickers]
R_sims_noarma <- R_sims_noarma[, tickers]

# ------------------------------------------------------------
# Helper to compute portfolio stats under MNTS sims
# ------------------------------------------------------------
portfolio_stats <- function(w, R_mat, alpha = 0.95) {
  port_ret <- as.numeric(R_mat %*% w)  # simulated 21d returns
  
  mean_ret <- mean(port_ret)
  sd_ret   <- stats::sd(port_ret)
  
  losses <- -port_ret
  var_a  <- unname(stats::quantile(losses, alpha, type = 7))
  cvar_a <- mean(losses[losses >= var_a])
  
  list(
    mean = mean_ret,
    sd   = sd_ret,
    var  = var_a,
    cvar = cvar_a,
    sims = port_ret
  )
}

# ------------------------------------------------------------
# 2) Load historical max-Sharpe weights
# ------------------------------------------------------------
w_ms_path <- "outputs/tables/weights_max_sharpe.csv"
if (!file.exists(w_ms_path)) {
  stop("Cannot find ", w_ms_path, " — run 04_max_sharpe.R first.")
}

w_ms_df <- read.csv(w_ms_path, stringsAsFactors = FALSE)

name_col <- if ("ticker" %in% names(w_ms_df)) {
  "ticker"
} else if ("symbol" %in% names(w_ms_df)) {
  "symbol"
} else {
  stop("weights_max_sharpe.csv must have 'ticker' or 'symbol' column.")
}

weight_col <- if ("weight" %in% names(w_ms_df)) {
  "weight"
} else if ("w" %in% names(w_ms_df)) {
  "w"
} else {
  stop("weights_max_sharpe.csv must have 'weight' or 'w' column.")
}

w_ms_vec <- setNames(w_ms_df[[weight_col]], w_ms_df[[name_col]])

missing_in_ms <- setdiff(tickers, names(w_ms_vec))
if (length(missing_in_ms) > 0) {
  stop("These tickers are in sims but not in weights_max_sharpe: ",
       paste(missing_in_ms, collapse = ", "))
}
w_ms_vec <- w_ms_vec[tickers]

# ------------------------------------------------------------
# 3) Load CVaR-frontier max-Sharpe weights (ARMA & no-ARMA)
# ------------------------------------------------------------
frontier_arma_path   <- "outputs/tables/frontier_cvar_21d.csv"
frontier_noarma_path <- "outputs/tables/frontier_cvar_21d_noarma.csv"

if (!file.exists(frontier_arma_path) || !file.exists(frontier_noarma_path)) {
  stop("Missing frontier_cvar csvs — run 07_frontier_cvar.R and 07b_frontier_cvar_noarma.R.")
}

frontier_arma   <- read.csv(frontier_arma_path,   stringsAsFactors = FALSE)
frontier_noarma <- read.csv(frontier_noarma_path, stringsAsFactors = FALSE)

if (!("Sharpe_sd" %in% names(frontier_arma)) ||
    !("Sharpe_sd" %in% names(frontier_noarma))) {
  stop("Both frontier tables must have a 'Sharpe_sd' column.")
}

best_idx_arma   <- which.max(frontier_arma$Sharpe_sd)
best_idx_noarma <- which.max(frontier_noarma$Sharpe_sd)

w_cols_arma   <- grep("^w_", names(frontier_arma),   value = TRUE)
w_cols_noarma <- grep("^w_", names(frontier_noarma), value = TRUE)

if (length(w_cols_arma) != n_assets || length(w_cols_noarma) != n_assets) {
  stop("Number of w_* columns in frontier tables does not match sims n_assets.")
}

w_cvar_arma   <- as.numeric(frontier_arma[best_idx_arma,   w_cols_arma])
w_cvar_noarma <- as.numeric(frontier_noarma[best_idx_noarma, w_cols_noarma])

names(w_cvar_arma)   <- sub("^w_", "", w_cols_arma)
names(w_cvar_noarma) <- sub("^w_", "", w_cols_noarma)

w_cvar_arma   <- w_cvar_arma[tickers]
w_cvar_noarma <- w_cvar_noarma[tickers]

# ------------------------------------------------------------
# 4) Predicted stats under MNTS sims (21d horizon)
# ------------------------------------------------------------
stats_ms          <- portfolio_stats(w_ms_vec,        R_sims_arma,   alpha)
stats_cvar_arma   <- portfolio_stats(w_cvar_arma,     R_sims_arma,   alpha)
stats_cvar_noarma <- portfolio_stats(w_cvar_noarma,   R_sims_noarma, alpha)

# ------------------------------------------------------------
# 5) Fetch REAL 2025 prices and build 21d windows
# ------------------------------------------------------------
last_train_date <- as.Date(cfg$end_date)          # should be 2024-12-31
backtest_start  <- as.Date("2025-01-01")
backtest_end    <- as.Date("2025-12-31")

if (backtest_start <= last_train_date) {
  stop("Config end_date must be <= 2024-12-31 to backtest full 2025.")
}

message("Downloading 2025 prices…")

prices_2025 <- tq_get(
  tickers,
  from = backtest_start - 15,   # buffer
  to   = backtest_end,
  get  = "stock.prices",
  complete_cases = FALSE
) %>%
  select(symbol, date, adjusted)

if (nrow(prices_2025) == 0) {
  stop("No price data returned for 2025. Check tickers or internet.")
}

daily_all <- prices_2025 %>%
  arrange(symbol, date) %>%
  group_by(symbol) %>%
  mutate(ret_log = log(adjusted / lag(adjusted))) %>%
  ungroup()

future_dates <- sort(unique(
  daily_all$date[daily_all$date >= backtest_start & daily_all$date <= backtest_end]
))

H        <- 21L   # 21 trading-day horizon
n_future <- length(future_dates)
n_windows <- n_future %/% H

if (n_windows < 1L) {
  stop("Not enough trading days in 2025 to form a single 21d window.")
}

message("Backtesting 2025 using ", n_windows,
        " non-overlapping 21d windows (", n_future, " trading days).")

real_ms_vec          <- numeric(n_windows)
real_cvar_arma_vec   <- numeric(n_windows)
real_cvar_noarma_vec <- numeric(n_windows)

window_start <- as.Date(rep(NA, n_windows))
window_end   <- as.Date(rep(NA, n_windows))

for (j in seq_len(n_windows)) {
  idx_start    <- (j - 1L) * H + 1L
  idx_end      <- idx_start + H - 1L
  window_dates <- future_dates[idx_start:idx_end]
  
  window_start[j] <- min(window_dates)
  window_end[j]   <- max(window_dates)
  
  rets_win <- daily_all %>%
    filter(date %in% window_dates) %>%
    select(date, symbol, ret_log) %>%
    pivot_wider(names_from = symbol, values_from = ret_log) %>%
    arrange(date) %>%
    drop_na()
  
  missing_cols <- setdiff(tickers, names(rets_win))
  if (length(missing_cols) > 0) {
    stop("Missing returns for some tickers in window ", j, ": ",
         paste(missing_cols, collapse = ", "))
  }
  
  R_daily_win <- as.matrix(rets_win[, tickers])
  
  ret_ms_daily          <- as.numeric(R_daily_win %*% w_ms_vec)
  ret_cvar_arma_daily   <- as.numeric(R_daily_win %*% w_cvar_arma)
  ret_cvar_noarma_daily <- as.numeric(R_daily_win %*% w_cvar_noarma)
  
  real_ms_vec[j]          <- sum(ret_ms_daily)
  real_cvar_arma_vec[j]   <- sum(ret_cvar_arma_daily)
  real_cvar_noarma_vec[j] <- sum(ret_cvar_noarma_daily)
}

# ------------------------------------------------------------
# 6) Percentiles & VaR breaches across windows
# ------------------------------------------------------------
sim_ms_21d          <- stats_ms$sims
sim_cvar_arma_21d   <- stats_cvar_arma$sims
sim_cvar_noarma_21d <- stats_cvar_noarma$sims

real_ms_loss_vec          <- -real_ms_vec
real_cvar_arma_loss_vec   <- -real_cvar_arma_vec
real_cvar_noarma_loss_vec <- -real_cvar_noarma_vec

# Percentile of realized 21d return within simulated distribution
perc_ms          <- sapply(real_ms_vec,          function(x) mean(sim_ms_21d          <= x))
perc_cvar_arma   <- sapply(real_cvar_arma_vec,   function(x) mean(sim_cvar_arma_21d   <= x))
perc_cvar_noarma <- sapply(real_cvar_noarma_vec, function(x) mean(sim_cvar_noarma_21d <= x))

# VaR(95) breach flags per window
VaR95_breach_ms          <- real_ms_loss_vec          > stats_ms$var
VaR95_breach_cvar_arma   <- real_cvar_arma_loss_vec   > stats_cvar_arma$var
VaR95_breach_cvar_noarma <- real_cvar_noarma_loss_vec > stats_cvar_noarma$var

# Detail per window
windows_df <- data.frame(
  window_id                = seq_len(n_windows),
  start_date               = window_start,
  end_date                 = window_end,
  realized_ms_21d          = real_ms_vec,
  realized_cvar_arma_21d   = real_cvar_arma_vec,
  realized_cvar_noarma_21d = real_cvar_noarma_vec,
  percentile_ms            = perc_ms,
  percentile_cvar_arma     = perc_cvar_arma,
  percentile_cvar_noarma   = perc_cvar_noarma,
  VaR95_breach_ms          = VaR95_breach_ms,
  VaR95_breach_cvar_arma   = VaR95_breach_cvar_arma,
  VaR95_breach_cvar_noarma = VaR95_breach_cvar_noarma,
  stringsAsFactors = FALSE
)

dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(
  windows_df,
  "outputs/tables/backtest_2025_windows_detail.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 7) Summarize + write CSV
# ------------------------------------------------------------
summary_df <- data.frame(
  portfolio              = c("MS_historical", "MS_CVaR_ARMA", "MS_CVaR_noARMA"),
  pred_mean_21d          = c(stats_ms$mean,
                             stats_cvar_arma$mean,
                             stats_cvar_noarma$mean),
  pred_sd_21d            = c(stats_ms$sd,
                             stats_cvar_arma$sd,
                             stats_cvar_noarma$sd),
  pred_VaR95_loss_21d    = c(stats_ms$var,
                             stats_cvar_arma$var,
                             stats_cvar_noarma$var),
  pred_CVaR95_loss_21d   = c(stats_ms$cvar,
                             stats_cvar_arma$cvar,
                             stats_cvar_noarma$cvar),
  realized_mean_21d      = c(mean(real_ms_vec),
                             mean(real_cvar_arma_vec),
                             mean(real_cvar_noarma_vec)),
  realized_sd_21d        = c(stats::sd(real_ms_vec),
                             stats::sd(real_cvar_arma_vec),
                             stats::sd(real_cvar_noarma_vec)),
  realized_mean_percentile = c(mean(perc_ms),
                               mean(perc_cvar_arma),
                               mean(perc_cvar_noarma)),
  n_windows              = n_windows,
  VaR95_breach_rate      = c(mean(VaR95_breach_ms),
                             mean(VaR95_breach_cvar_arma),
                             mean(VaR95_breach_cvar_noarma)),
  stringsAsFactors = FALSE
)

out_path <- "outputs/tables/backtest_2025_summary.csv"
write.csv(summary_df, out_path, row.names = FALSE)

message("Wrote summary to: ", out_path)
message("=== Done: 09_compare_2025.R ===")
