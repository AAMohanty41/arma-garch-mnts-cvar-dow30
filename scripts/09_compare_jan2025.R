#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(tidyquant)
  library(dplyr)
  library(tidyr)
  library(readr)
})

message("=== 09_compare_jan2025.R : compare predicted vs realized Jan 2025 ===")

# ------------------------------------------------------------
# 1) Load config + MNTS sims
# ------------------------------------------------------------

cfg <- yaml::read_yaml("config/config.yml")
tickers <- cfg$universe

sim_path <- "data/processed/sims_mnts_21d.rds"
if (!file.exists(sim_path)) {
  stop("Cannot find ", sim_path, " — run 06_mnts_fit_sim.R first.")
}

sims <- readRDS(sim_path)
R_sims <- as.matrix(sims)

S        <- nrow(R_sims)
n_assets <- ncol(R_sims)
sim_tickers  <- colnames(R_sims)

if (is.null(sim_tickers)) {
  stop("sims_mnts_21d.rds has no column names; they must match tickers.")
}

if (!setequal(sim_tickers, tickers)) {
  warning("Tickers in sims and config$universe differ; using intersection.")
}
tickers <- sim_tickers  # enforce same order as sims

alpha <- 0.95

portfolio_stats <- function(w, R_mat, alpha = 0.95) {
  port_ret <- as.numeric(R_mat %*% w)   # 21d simulated returns
  
  mean_ret <- mean(port_ret)
  sd_ret   <- stats::sd(port_ret)
  
  losses <- -port_ret
  var_alpha  <- unname(stats::quantile(losses, alpha, type = 7))
  cvar_alpha <- mean(losses[losses >= var_alpha])
  
  list(
    mean = mean_ret,
    sd   = sd_ret,
    var  = var_alpha,
    cvar = cvar_alpha,
    sims = port_ret       # for percentiles
  )
}

# ------------------------------------------------------------
# 2) Historical max-Sharpe weights
# ------------------------------------------------------------

w_ms_path <- "outputs/tables/weights_max_sharpe.csv"
if (!file.exists(w_ms_path)) {
  stop("Cannot find ", w_ms_path, " — run 04_max_sharpe.R first.")
}





w_ms_df <- read.csv(w_ms_path, stringsAsFactors = FALSE)

# Figure out which column is the name and which is the weight
name_col <- NULL
if ("ticker" %in% names(w_ms_df)) {
  name_col <- "ticker"
} else if ("symbol" %in% names(w_ms_df)) {
  name_col <- "symbol"
} else {
  stop("weights_max_sharpe.csv must have a 'ticker' or 'symbol' column; found: ",
       paste(names(w_ms_df), collapse = ", "))
}

weight_col <- NULL
if ("weight" %in% names(w_ms_df)) {
  weight_col <- "weight"
} else if ("w" %in% names(w_ms_df)) {
  weight_col <- "w"
} else {
  stop("weights_max_sharpe.csv must have a 'weight' or 'w' column; found: ",
       paste(names(w_ms_df), collapse = ", "))
}

w_ms_vec <- setNames(w_ms_df[[weight_col]], w_ms_df[[name_col]])





missing_in_ms <- setdiff(tickers, names(w_ms_vec))
if (length(missing_in_ms) > 0) {
  stop("These tickers are in sims but not in weights_max_sharpe: ",
       paste(missing_in_ms, collapse = ", "))
}
w_ms_vec <- w_ms_vec[tickers]

# ------------------------------------------------------------
# 3) CVaR-frontier max-Sharpe weights
# ------------------------------------------------------------

frontier_path <- "outputs/tables/frontier_cvar_21d.csv"
if (!file.exists(frontier_path)) {
  stop("Cannot find ", frontier_path, " — run 07_frontier_cvar.R first.")
}

frontier <- read.csv(frontier_path, stringsAsFactors = FALSE)

if (!("Sharpe_sd" %in% names(frontier))) {
  stop("frontier_cvar_21d.csv has no 'Sharpe_sd' column.")
}

best_idx <- which.max(frontier$Sharpe_sd)

w_cols <- grep("^w_", names(frontier), value = TRUE)
if (length(w_cols) != n_assets) {
  stop("Number of w_* columns in frontier_cvar_21d.csv (", length(w_cols),
       ") does not match assets in sims (", n_assets, ").")
}

w_cvar_vec <- as.numeric(frontier[best_idx, w_cols])
names(w_cvar_vec) <- sub("^w_", "", w_cols)
w_cvar_vec <- w_cvar_vec[tickers]

# ------------------------------------------------------------
# 4) Predicted stats under MNTS sims
# ------------------------------------------------------------

stats_ms   <- portfolio_stats(w_ms_vec,   R_sims, alpha)
stats_cvar <- portfolio_stats(w_cvar_vec, R_sims, alpha)




# ------------------------------------------------------------
# 5) Fetch REAL early-2025 prices and compute 21-day returns
#    Use the first 21 trading days AFTER the training window.
# ------------------------------------------------------------

message("Downloading early-2025 prices…")

last_train_date <- as.Date(cfg$end_date)  # e.g. "2024-12-31"

prices_early <- tq_get(
  tickers,
  from = last_train_date - 15,   # a bit before, for lag
  to   = last_train_date + 60,   # ~3 months after, enough buffer
  get  = "stock.prices",
  complete_cases = FALSE
) %>%
  select(symbol, date, adjusted)

if (nrow(prices_early) == 0) {
  stop("No price data returned for early 2025. Check internet / ticker availability.")
}

daily_all <- prices_early %>%
  arrange(symbol, date) %>%
  group_by(symbol) %>%
  mutate(
    ret_log = log(adjusted / lag(adjusted))
  ) %>%
  ungroup()

# First 21 trading days AFTER the training window
future_dates <- sort(unique(daily_all$date[daily_all$date > last_train_date]))
h <- 21
if (length(future_dates) < h) {
  stop("Fewer than ", h, " trading days after ", last_train_date,
       " in fetched data. Got: ", length(future_dates))
}

window_dates <- future_dates[1:h]

rets_win <- daily_all %>%
  filter(date %in% window_dates) %>%
  select(date, symbol, ret_log) %>%
  pivot_wider(names_from = symbol, values_from = ret_log) %>%
  arrange(date) %>%
  drop_na()


missing_cols <- setdiff(tickers, names(rets_win))
if (length(missing_cols) > 0) {
  stop("These tickers are in sims but missing returns in Jan 2025: ",
       paste(missing_cols, collapse = ", "))
}

R_daily <- as.matrix(rets_win[, tickers])

# Portfolio daily returns and 21-day totals
ret_ms_daily   <- as.numeric(R_daily %*% w_ms_vec)
ret_cvar_daily <- as.numeric(R_daily %*% w_cvar_vec)

real_ms_21d   <- sum(ret_ms_daily)
real_cvar_21d <- sum(ret_cvar_daily)



# ------------------------------------------------------------
# 6) Percentiles & VaR breaches
# ------------------------------------------------------------

sim_ms_21d   <- stats_ms$sims
sim_cvar_21d <- stats_cvar$sims

perc_ms   <- mean(sim_ms_21d   <= real_ms_21d)
perc_cvar <- mean(sim_cvar_21d <= real_cvar_21d)

real_ms_loss   <- -real_ms_21d
real_cvar_loss <- -real_cvar_21d

summary_df <- data.frame(
  portfolio           = c("MS_historical", "MS_CVaR_frontier"),
  pred_mean_21d       = c(stats_ms$mean,   stats_cvar$mean),
  pred_sd_21d         = c(stats_ms$sd,     stats_cvar$sd),
  pred_VaR95_loss     = c(stats_ms$var,    stats_cvar$var),
  pred_CVaR95_loss    = c(stats_ms$cvar,   stats_cvar$cvar),
  realized_21d        = c(real_ms_21d,     real_cvar_21d),
  realized_loss       = c(real_ms_loss,    real_cvar_loss),
  realized_percentile = c(perc_ms,         perc_cvar),
  VaR95_breach        = c(real_ms_loss   > stats_ms$var,
                          real_cvar_loss > stats_cvar$var),
  stringsAsFactors = FALSE
)
jan2025_backtest_summary <- summary_df   # put in Environment when sourced

dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
out_path <- "outputs/tables/jan2025_backtest_summary.csv"
write.csv(summary_df, out_path, row.names = FALSE)

message("Wrote summary to: ", out_path)
message("=== Done: 09_compare_jan2025.R ===")
